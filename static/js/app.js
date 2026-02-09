/**
 * app.js - Loads data from static JSON files (GitHub Pages / CDN compatible)
 *
 * For GitHub Pages: DATA_BASE = './data'
 * For Bunny CDN:    DATA_BASE = 'https://your-zone.b-cdn.net'
 */

// ============= CONFIGURATION =============
// Change this URL when moving to Bunny CDN
const DATA_BASE = './data';
// =========================================

// Get pair ID from URL
const urlParams = new URLSearchParams(window.location.search);
const PAIR_ID = urlParams.get('pair');

// Family index cache (loaded once)
let FAMILY_INDEX = null;

// These will be populated from API
let DATA = null;
let SUMMARY = null;
let PDB64_FULL = "";
let PLMA_DATA = null;

// PDBe highlight state
let isProteinHighlighted = false;

// Block Molstar volume server
(function() {
    const block = ['molstarvolseg.ncbr.muni.cz', 'localhost:9000'];
    const _f = window.fetch.bind(window);
    window.fetch = (u, i) => {
        try { if (block.some(h => (typeof u === 'string' ? u : u?.url || '').includes(h))) return Promise.resolve(new Response('{"items":[]}', {status:200})); } catch(e){}
        return _f(u, i);
    };
})();

// Main initialization - loads data then runs notebook code
async function loadDataAndInit() {
    if (!PAIR_ID) {
        document.getElementById('loadingOverlay').textContent = 'No pair specified';
        return;
    }
    
    try {
        // Static file fetches (GitHub Pages / CDN compatible)
        const [dataResp, summaryResp, variantsResp] = await Promise.all([
            fetch(`${DATA_BASE}/pairs/${PAIR_ID}/report.json`),
            fetch(`${DATA_BASE}/pairs/${PAIR_ID}/summary.json`),
            fetch(`${DATA_BASE}/pairs/${PAIR_ID}/pdb.json`)
        ]);

        if (!dataResp.ok) throw new Error(`Failed to load pair: ${PAIR_ID}`);

        DATA = await dataResp.json();
        SUMMARY = summaryResp.ok ? await summaryResp.json() : { gene1: {}, gene2: {}, pair: {}, conservation: {}, boxplots: {} };

        // Load PDB variants from combined pdb.json file
        if (variantsResp.ok) {
            const variants = await variantsResp.json();
            PDB64_FULL = variants.pdb64_full || "";
            window.PDB64_A = variants.pdb64_a || PDB64_FULL;
            window.PDB64_B = variants.pdb64_b || PDB64_FULL;

            // Store AM-colored PDB variants by mode
            window.PDB64_AM_BY_MODE = {};
            for (const key in variants) {
                if (key.startsWith('pdb64_am_')) {
                    const mode = key.replace('pdb64_am_', '');
                    window.PDB64_AM_BY_MODE[mode] = variants[key];
                }
            }

            // Store other color mode variants
            window.PDB64_PLDDT = variants.pdb64_plddt || PDB64_FULL;
            window.PDB64_ALIGNED = variants.pdb64_aligned || null;
            window.PDB64_DOMAINS = variants.pdb64_domains || null;

            console.log('Loaded PDB variants:', Object.keys(variants));
        } else {
            window.PDB64_A = PDB64_FULL;
            window.PDB64_B = PDB64_FULL;
            window.PDB64_AM_BY_MODE = {};
            window.PDB64_PLDDT = PDB64_FULL;
            window.PDB64_ALIGNED = null;
            window.PDB64_DOMAINS = null;
        }
        
        // Initialize DATA-dependent variables
        AM_MODES = DATA.amModes || ['raw'];
        PDBe_COMPLEXES = DATA.pdbeComplexes || [];
        UNIPROT_A = DATA.a1 || '';
        UNIPROT_B = DATA.a2 || '';
        
        // Update title
        document.title = `${DATA.g1} vs ${DATA.g2}`;
        document.getElementById('titleMain').textContent = `${DATA.g1} ↔ ${DATA.g2}`;
        document.getElementById('titleSub').textContent = `Paralog pair ${DATA.PAIR}`;

        // Load PLMA alignment data
        try {
            const plmaResp = await fetch(`${DATA_BASE}/pairs/${PAIR_ID}/plma.json`);
            if (plmaResp.ok) PLMA_DATA = await plmaResp.json();
        } catch(e) { console.warn('PLMA data not available:', e); }

        // Load family data
        await loadFamilyData();

        // Now run the notebook initialization code (main is defined at the end of the file)
        await main();
        
        document.getElementById('loadingOverlay').style.display = 'none';
    } catch(e) {
        console.error(e);
        document.getElementById('loadingOverlay').textContent = 'Error: ' + e.message;
    }
}

// ========== FAMILY CONSTELLATION VIEW ==========
// State for the constellation
let constellationState = {
  centerGene: null,           // Gene at center of constellation
  selectedGenes: [],          // Currently selected genes (max 2)
  hoveredGene: null,          // Gene being hovered
  allGenes: [],               // All genes in the family
  geneData: {},               // Gene -> { pairs, identities, hasData }
  pairData: {},               // pairId -> pair metadata
  closestPair: null,          // The pair with highest identity
};

// Full families data cache
let FULL_FAMILIES = null;

// Family navigation
async function loadFamilyData() {
  if (!DATA || (!DATA.g1 && !DATA.g2)) {
    console.log('Family: No DATA or genes');
    return;
  }

  const familyNav = document.getElementById('familyNav');
  if (!familyNav) {
    console.log('Family: DOM elements not found');
    return;
  }

  try {
    console.log(`Loading family data for ${DATA.g1} and ${DATA.g2}...`);

    // Load full families data (includes all genes in each family with identities)
    if (!FULL_FAMILIES) {
      const fullFamiliesResp = await fetch(`${DATA_BASE}/full_families.json`);
      if (fullFamiliesResp.ok) {
        FULL_FAMILIES = await fullFamiliesResp.json();
        console.log('Loaded full_families.json');
      }
    }

    // Load family index (for pairs with reports)
    if (!FAMILY_INDEX) {
      const indexResp = await fetch(`${DATA_BASE}/family_index.json`);
      if (indexResp.ok) {
        FAMILY_INDEX = await indexResp.json();
      }
    }

    // Load index.json to get all pair metadata (for pairs with reports)
    const indexResp = await fetch(`${DATA_BASE}/index.json`);
    const allPairsIndex = indexResp.ok ? await indexResp.json() : [];
    const pairMap = new Map(allPairsIndex.map(p => [p.id, p]));

    // Get full family members from full_families.json
    let familyGenes = new Set([DATA.g1, DATA.g2]);
    let fullFamilyIdentities = {};

    if (FULL_FAMILIES) {
      // Find which family these genes belong to
      const familyId = FULL_FAMILIES.families[DATA.g1] || FULL_FAMILIES.families[DATA.g2];
      if (familyId && FULL_FAMILIES.family_data[familyId]) {
        const familyData = FULL_FAMILIES.family_data[familyId];
        // Add all genes from the full family
        familyData.genes.forEach(g => familyGenes.add(g));
        fullFamilyIdentities = familyData.identities || {};
        console.log(`Full family ${familyId}: ${familyData.genes.length} genes`);
      }
    }

    // Also add genes from pairs with reports (in case full_families is missing some)
    const allPairIds = new Set();

    // Get pairs for g1
    (FAMILY_INDEX && FAMILY_INDEX[DATA.g1] || []).forEach(pairId => {
      allPairIds.add(pairId);
      const meta = pairMap.get(pairId);
      if (meta) {
        familyGenes.add(meta.geneA);
        familyGenes.add(meta.geneB);
      }
    });

    // Get pairs for g2
    (FAMILY_INDEX && FAMILY_INDEX[DATA.g2] || []).forEach(pairId => {
      allPairIds.add(pairId);
      const meta = pairMap.get(pairId);
      if (meta) {
        familyGenes.add(meta.geneA);
        familyGenes.add(meta.geneB);
      }
    });

    // Build gene data with identities to other genes
    const geneData = {};
    familyGenes.forEach(gene => {
      // Get identities from full_families.json (sequence identity from CSV)
      const fullIdentities = fullFamilyIdentities[gene] || {};

      geneData[gene] = {
        gene: gene,
        pairs: [],        // Pairs this gene is part of (with reports)
        identities: { ...fullIdentities },  // Start with full family identities
        hasData: false    // Has at least one pair with a report
      };
    });

    // Process pairs with reports to update hasData and override identities if needed
    const pairData = {};
    let closestPair = null;
    let maxIdent = -1;

    allPairIds.forEach(pairId => {
      const meta = pairMap.get(pairId);
      if (!meta) return;

      const pair = {
        pair_id: pairId,
        gene_a: meta.geneA,
        gene_b: meta.geneB,
        fident: meta.fident,
        tm_score: meta.tm,
        has_report: true
      };
      pairData[pairId] = pair;

      // Update gene data - mark as having data
      if (geneData[pair.gene_a]) {
        geneData[pair.gene_a].pairs.push(pairId);
        // Use structural identity (fident) from report if available
        if (pair.fident) {
          geneData[pair.gene_a].identities[pair.gene_b] = pair.fident;
        }
        geneData[pair.gene_a].hasData = true;
      }
      if (geneData[pair.gene_b]) {
        geneData[pair.gene_b].pairs.push(pairId);
        if (pair.fident) {
          geneData[pair.gene_b].identities[pair.gene_a] = pair.fident;
        }
        geneData[pair.gene_b].hasData = true;
      }

      // Track closest pair (only among those with reports)
      if (pair.fident && pair.fident > maxIdent) {
        maxIdent = pair.fident;
        closestPair = pair;
      }
    });

    // Also check full_families for the closest pair if we don't have one
    if (!closestPair && Object.keys(fullFamilyIdentities).length > 0) {
      // Find the pair with highest identity in the full family
      Object.entries(fullFamilyIdentities).forEach(([geneA, identities]) => {
        Object.entries(identities).forEach(([geneB, ident]) => {
          if (ident > maxIdent) {
            maxIdent = ident;
            closestPair = { gene_a: geneA, gene_b: geneB, fident: ident, has_report: false };
          }
        });
      });
    }

    // Update constellation state
    constellationState = {
      centerGene: DATA.g1,  // Start with g1 as center
      selectedGenes: [DATA.g1, DATA.g2],
      hoveredGene: null,
      allGenes: Array.from(familyGenes),
      geneData: geneData,
      pairData: pairData,
      closestPair: closestPair
    };

    const pairsWithReports = allPairIds.size;
    console.log(`Family: ${familyGenes.size} genes total, ${pairsWithReports} pairs with reports`);

    // Show the family section
    familyNav.style.display = 'block';
    const familySubtitle = document.getElementById('familySubtitle');
    if (familySubtitle) {
      familySubtitle.textContent = `${familyGenes.size} genes in family · ${pairsWithReports} pairs with reports`;
    }

    // Initialize constellation canvas
    initFamilyConstellation();

  } catch (e) {
    console.error('Failed to load family data:', e);
  }
}

function initFamilyConstellation() {
  const canvas = document.getElementById('familyNetworkCanvas');
  if (!canvas) return;

  // Handle HiDPI displays for sharp rendering - use minimum 2x for crisp rendering
  const dpr = Math.max(2, window.devicePixelRatio || 1);
  const container = canvas.parentElement;
  // Allow larger canvas on big screens for better resolution
  const displayWidth = Math.min(container.offsetWidth - 20, 1200);
  const displayHeight = Math.min(560, Math.max(400, displayWidth * 0.5));

  // Set canvas size accounting for device pixel ratio (minimum 2x for crisp text)
  canvas.width = displayWidth * dpr;
  canvas.height = displayHeight * dpr;
  canvas.style.width = displayWidth + 'px';
  canvas.style.height = displayHeight + 'px';

  // Scale context to match
  const ctx = canvas.getContext('2d');
  ctx.scale(dpr, dpr);

  // Store display dimensions for calculations
  canvas._displayWidth = displayWidth;
  canvas._displayHeight = displayHeight;
  canvas._dpr = dpr;

  // Add event listeners
  canvas.removeEventListener('mousemove', handleConstellationMouseMove);
  canvas.removeEventListener('click', handleConstellationClick);
  canvas.addEventListener('mousemove', handleConstellationMouseMove);
  canvas.addEventListener('click', handleConstellationClick);

  // Initial render
  renderFamilyConstellation();
}

function handleConstellationMouseMove(e) {
  const canvas = e.target;
  const rect = canvas.getBoundingClientRect();
  const x = e.clientX - rect.left;
  const y = e.clientY - rect.top;

  const hoveredGene = findGeneAtPosition(x, y);
  if (hoveredGene !== constellationState.hoveredGene) {
    constellationState.hoveredGene = hoveredGene;
    canvas.style.cursor = hoveredGene ? 'pointer' : 'default';
    renderFamilyConstellation();
  }
}

function handleConstellationClick(e) {
  const canvas = e.target;
  const rect = canvas.getBoundingClientRect();
  const x = e.clientX - rect.left;
  const y = e.clientY - rect.top;

  const clickedGene = findGeneAtPosition(x, y);
  if (!clickedGene) return;

  const state = constellationState;
  const geneInfo = state.geneData[clickedGene];

  // Check if clicked gene is currently selected (for unselect logic)
  const isCenterAndSelected = clickedGene === state.centerGene && state.selectedGenes.includes(clickedGene);
  const isPairPartner = state.selectedGenes.length === 2 && state.selectedGenes[1] === clickedGene;

  // CASE 1: Clicking a SELECTED center gene - unselect it
  if (isCenterAndSelected) {
    if (state.selectedGenes.length === 2) {
      // Had pair partner - partner becomes new center, constellation rebuilds
      const oldPartner = state.selectedGenes[1];
      state.centerGene = oldPartner;
      state.selectedGenes = [oldPartner];
      renderFamilyConstellation();
    } else {
      // Only center selected - clear selection but keep centerGene for positioning
      // Visual shows nothing selected, next click will set new center
      state.selectedGenes = [];
      renderFamilyConstellation();
    }
    return;
  }

  // CASE 2: Clicking the pair partner - unselect it
  if (isPairPartner) {
    state.selectedGenes = [state.centerGene];
    renderFamilyConstellation();
    return;
  }

  // CASE 3: Clicking an unselected gene - check hasData first
  if (!geneInfo?.hasData) {
    console.log(`${clickedGene} has no report data - cannot select`);
    return;
  }

  // No selection at all - make this the new center and rebuild constellation
  if (state.selectedGenes.length === 0) {
    state.centerGene = clickedGene;
    state.selectedGenes = [clickedGene];
    renderFamilyConstellation();
    return;
  }

  // Only center selected - add as pair partner and load pair
  if (state.selectedGenes.length === 1) {
    state.selectedGenes = [state.centerGene, clickedGene];
    const pairId = findPairBetweenGenes(state.centerGene, clickedGene);
    if (pairId && state.pairData[pairId]) {
      window.location.href = `report.html?pair=${pairId}`;
    } else {
      renderFamilyConstellation();
    }
    return;
  }

  // Both selected - replace pair partner with new gene and load pair
  state.selectedGenes = [state.centerGene, clickedGene];
  const pairId = findPairBetweenGenes(state.centerGene, clickedGene);
  if (pairId && state.pairData[pairId]) {
    window.location.href = `report.html?pair=${pairId}`;
  } else {
    renderFamilyConstellation();
  }
}

function findPairBetweenGenes(gene1, gene2) {
  const pairId1 = `${gene1}_${gene2}`;
  const pairId2 = `${gene2}_${gene1}`;
  if (constellationState.pairData[pairId1]) return pairId1;
  if (constellationState.pairData[pairId2]) return pairId2;
  return null;
}

function findGeneAtPosition(x, y) {
  const positions = calculateGenePositions();
  const hitRadius = 20;

  for (const [gene, pos] of Object.entries(positions)) {
    const dx = x - pos.x;
    const dy = y - pos.y;
    if (dx * dx + dy * dy < hitRadius * hitRadius) {
      return gene;
    }
  }
  return null;
}

function getCanvasDisplayDimensions(canvas) {
  // Get the actual display dimensions (CSS pixels, not device pixels)
  // Priority: stored _displayWidth, then style.width, then getBoundingClientRect
  if (canvas._displayWidth && canvas._displayHeight) {
    return { width: canvas._displayWidth, height: canvas._displayHeight };
  }
  if (canvas.style.width && canvas.style.height) {
    return {
      width: parseInt(canvas.style.width, 10),
      height: parseInt(canvas.style.height, 10)
    };
  }
  const rect = canvas.getBoundingClientRect();
  return { width: rect.width, height: rect.height };
}

function calculateGenePositions() {
  const canvas = document.getElementById('familyNetworkCanvas');
  if (!canvas) return {};

  const state = constellationState;
  // Use display dimensions for calculations (CSS pixels, not device pixels)
  const { width, height } = getCanvasDisplayDimensions(canvas);
  const cx = width / 2;
  const cy = height / 2;
  // Leave generous padding for labels
  const maxRadius = Math.min(width, height) / 2 - 60;

  const positions = {};
  const centerGene = state.centerGene;

  // Center gene position
  if (centerGene) {
    positions[centerGene] = { x: cx, y: cy };
  }

  // Get other genes and compute their distance to center
  const otherGenes = state.allGenes.filter(g => g !== centerGene);
  if (otherGenes.length === 0) return positions;

  // Get identity to center for each gene
  const geneDistances = otherGenes.map(gene => {
    const identity = state.geneData[centerGene]?.identities[gene] || 0;
    return { gene, identity, distanceToCenter: 1 - identity };
  });

  // Angular clustering: order genes so similar genes are adjacent
  // Build inter-gene similarity lookup
  const getInterGeneIdentity = (g1, g2) => {
    return state.geneData[g1]?.identities[g2] ||
           state.geneData[g2]?.identities[g1] || 0;
  };

  // Greedy ordering: start with gene closest to center, then add most similar unplaced gene
  const ordered = [];
  const remaining = new Set(otherGenes);

  // Start with gene closest to center (highest identity)
  geneDistances.sort((a, b) => a.distanceToCenter - b.distanceToCenter);
  const firstGene = geneDistances[0].gene;
  ordered.push(firstGene);
  remaining.delete(firstGene);

  // Greedily add genes by similarity to the last placed gene
  while (remaining.size > 0) {
    const lastPlaced = ordered[ordered.length - 1];
    let bestGene = null;
    let bestSimilarity = -1;

    for (const gene of remaining) {
      const sim = getInterGeneIdentity(lastPlaced, gene);
      if (sim > bestSimilarity) {
        bestSimilarity = sim;
        bestGene = gene;
      }
    }

    // If no similarity found, just pick the one closest to center
    if (bestGene === null || bestSimilarity === 0) {
      let minDist = Infinity;
      for (const gene of remaining) {
        const dist = 1 - (state.geneData[centerGene]?.identities[gene] || 0);
        if (dist < minDist) {
          minDist = dist;
          bestGene = gene;
        }
      }
    }

    if (bestGene) {
      ordered.push(bestGene);
      remaining.delete(bestGene);
    } else {
      break;
    }
  }

  // Position genes in a circle, with radial distance based on identity to center
  // Angular position based on clustering order
  ordered.forEach((gene, i) => {
    const identity = state.geneData[centerGene]?.identities[gene] || 0;
    // Radial distance: higher identity = closer to center
    // Clamp minimum radius so genes don't overlap with center
    const minRadius = 40;
    const radius = minRadius + (1 - identity) * (maxRadius - minRadius);

    // Angular position based on order
    const angle = (i / ordered.length) * 2 * Math.PI - Math.PI / 2;

    positions[gene] = {
      x: cx + radius * Math.cos(angle),
      y: cy + radius * Math.sin(angle)
    };
  });

  return positions;
}

function renderFamilyConstellation() {
  const canvas = document.getElementById('familyNetworkCanvas');
  if (!canvas) return;

  const ctx = canvas.getContext('2d');
  const dpr = canvas._dpr || window.devicePixelRatio || 1;

  // Use display dimensions (CSS pixels, not device pixels)
  const { width, height } = getCanvasDisplayDimensions(canvas);
  const cx = width / 2;
  const cy = height / 2;
  // Match the maxRadius calculation in calculateGenePositions
  const maxRadius = Math.min(width, height) / 2 - 60;

  const state = constellationState;
  const positions = calculateGenePositions();

  // Reset transform and clear canvas
  ctx.setTransform(1, 0, 0, 1, 0, 0);
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  // Apply DPR scaling for all subsequent drawing
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

  // Light background matching rest of page
  ctx.fillStyle = '#fafafa';
  ctx.fillRect(0, 0, width, height);

  // Draw orbit circles (identity thresholds) - including 0%
  const thresholds = [0, 0.2, 0.4, 0.6, 0.8];
  const minRadius = 40;
  thresholds.forEach(threshold => {
    const radius = minRadius + (1 - threshold) * (maxRadius - minRadius);
    ctx.beginPath();
    ctx.arc(cx, cy, radius, 0, 2 * Math.PI);
    ctx.strokeStyle = threshold === 0 ? 'rgba(0,0,0,0.2)' : 'rgba(0,0,0,0.1)';
    ctx.lineWidth = threshold === 0 ? 1.5 : 1;
    ctx.setLineDash([4, 6]);
    ctx.stroke();
    ctx.setLineDash([]);

    // Label the orbit
    ctx.font = '10px -apple-system, sans-serif';
    ctx.textAlign = 'left';
    ctx.fillStyle = '#666';
    ctx.fillText(`${(threshold * 100).toFixed(0)}%`, cx + radius + 5, cy - 2);
  });

  // Draw center crosshair
  ctx.beginPath();
  ctx.moveTo(cx - 6, cy);
  ctx.lineTo(cx + 6, cy);
  ctx.moveTo(cx, cy - 6);
  ctx.lineTo(cx, cy + 6);
  ctx.strokeStyle = 'rgba(0,0,0,0.15)';
  ctx.lineWidth = 1;
  ctx.stroke();

  // Find closest gene to center (highest identity)
  let closestToCenter = null;
  let highestIdentity = 0;
  if (state.centerGene && state.geneData[state.centerGene]) {
    const centerIdentities = state.geneData[state.centerGene].identities || {};
    for (const [gene, identity] of Object.entries(centerIdentities)) {
      if (identity > highestIdentity) {
        highestIdentity = identity;
        closestToCenter = gene;
      }
    }
  }

  // Draw double edge to closest gene (always visible)
  if (state.centerGene && closestToCenter && positions[closestToCenter] && positions[state.centerGene]) {
    const p1 = positions[state.centerGene];
    const p2 = positions[closestToCenter];
    const dx = p2.x - p1.x;
    const dy = p2.y - p1.y;
    const len = Math.sqrt(dx * dx + dy * dy);
    if (len > 0) {
      const nx = -dy / len * 3; // perpendicular offset
      const ny = dx / len * 3;

      // Check if closest pair matches selected pair
      const isSelectedPair = state.selectedGenes.length === 2 &&
        state.selectedGenes.includes(state.centerGene) &&
        state.selectedGenes.includes(closestToCenter);

      ctx.beginPath();
      ctx.moveTo(p1.x + nx, p1.y + ny);
      ctx.lineTo(p2.x + nx, p2.y + ny);
      ctx.moveTo(p1.x - nx, p1.y - ny);
      ctx.lineTo(p2.x - nx, p2.y - ny);
      ctx.strokeStyle = isSelectedPair ? '#0d9488' : '#888';
      ctx.lineWidth = 2;
      ctx.stroke();
    }
  }

  // Draw edge for selected pair (teal line) - only when BOTH genes are selected
  // Skip if it's already the closest pair (already drawn above)
  const bothSelected = state.selectedGenes.length === 2;
  const selectedIsClosest = bothSelected &&
    state.selectedGenes.includes(state.centerGene) &&
    state.selectedGenes.includes(closestToCenter);
  if (bothSelected && !selectedIsClosest && positions[state.selectedGenes[0]] && positions[state.selectedGenes[1]]) {
    const p1 = positions[state.selectedGenes[0]];
    const p2 = positions[state.selectedGenes[1]];
    ctx.beginPath();
    ctx.moveTo(p1.x, p1.y);
    ctx.lineTo(p2.x, p2.y);
    ctx.strokeStyle = '#0d9488';
    ctx.lineWidth = 3;
    ctx.stroke();
  }

  // Draw genes (nodes)
  state.allGenes.forEach(gene => {
    const pos = positions[gene];
    if (!pos) return;

    const geneInfo = state.geneData[gene];
    const isPositionalCenter = gene === state.centerGene;
    const isInSelection = state.selectedGenes.includes(gene);
    // Visual "center" only if both positional center AND in selection
    const isSelectedCenter = isPositionalCenter && isInSelection;
    // Partner is the second gene in selection (not the center)
    const isSelectedPartner = state.selectedGenes.length === 2 && state.selectedGenes[1] === gene;
    const isHovered = gene === state.hoveredGene;
    const hasData = geneInfo?.hasData;

    // Determine colors and sizes - scheme for light background
    let fillColor, strokeColor, radius, labelColor;

    if (isSelectedCenter) {
      // Selected center gene - warm amber/orange
      fillColor = '#d97706';
      strokeColor = '#92400e';
      radius = 14;
      labelColor = '#92400e';
    } else if (isSelectedPartner) {
      // Selected partner gene - purple
      fillColor = '#7c3aed';
      strokeColor = '#5b21b6';
      radius = 12;
      labelColor = '#5b21b6';
    } else if (hasData) {
      // Genes with data - slate blue
      fillColor = isHovered ? '#475569' : '#64748b';
      strokeColor = isHovered ? '#334155' : '#475569';
      radius = isHovered ? 10 : 8;
      labelColor = '#334155';
    } else {
      // No data - light gray
      fillColor = '#d1d5db';
      strokeColor = '#9ca3af';
      radius = 7;
      labelColor = '#9ca3af';
    }

    // Draw node glow for selected genes
    if (isSelectedCenter || isSelectedPartner) {
      ctx.beginPath();
      ctx.arc(pos.x, pos.y, radius + 6, 0, 2 * Math.PI);
      ctx.fillStyle = fillColor + '30';
      ctx.fill();
    }

    // Draw node
    ctx.beginPath();
    ctx.arc(pos.x, pos.y, radius, 0, 2 * Math.PI);
    ctx.fillStyle = fillColor;
    ctx.fill();
    ctx.strokeStyle = strokeColor;
    ctx.lineWidth = 2;
    ctx.stroke();

    // Draw label
    ctx.fillStyle = labelColor;
    ctx.font = (isSelectedCenter || isSelectedPartner) ? 'bold 11px sans-serif' : '10px sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';

    // Position label below node
    const labelY = pos.y + radius + 12;
    ctx.fillText(gene, pos.x, labelY);

    // Show identity on hover (not for gene at center position)
    if (isHovered && !isPositionalCenter && state.centerGene) {
      const identity = state.geneData[state.centerGene]?.identities[gene];
      if (identity !== undefined) {
        ctx.fillStyle = '#d97706';
        ctx.font = 'bold 10px sans-serif';
        ctx.fillText(`${(identity * 100).toFixed(1)}%`, pos.x, labelY + 11);
      }
    }
  });

  // Draw legend
  drawConstellationLegend(ctx, width, height);
}

function drawConstellationLegend(ctx, width, height) {
  const legendX = 15;
  const legendY = height - 130;

  ctx.font = '10px sans-serif';
  ctx.textAlign = 'left';
  ctx.textBaseline = 'middle';

  const items = [
    { color: '#d97706', label: 'Center gene (query)' },
    { color: '#0d9488', label: 'Current pair' },
    { color: '#64748b', label: 'Family member' },
    { color: '#d1d5db', label: 'No report data' },
  ];

  items.forEach((item, i) => {
    const y = legendY + i * 16;

    ctx.beginPath();
    ctx.arc(legendX + 5, y, 4, 0, 2 * Math.PI);
    ctx.fillStyle = item.color;
    ctx.fill();
    ctx.strokeStyle = '#666';
    ctx.lineWidth = 1;
    ctx.stroke();

    ctx.fillStyle = '#444';
    ctx.fillText(item.label, legendX + 16, y);
  });

  // Edge legend - current pair (teal)
  ctx.beginPath();
  ctx.moveTo(legendX, legendY + 70);
  ctx.lineTo(legendX + 12, legendY + 70);
  ctx.strokeStyle = '#0d9488';
  ctx.lineWidth = 3;
  ctx.stroke();

  ctx.fillStyle = '#444';
  ctx.fillText('Shown pair', legendX + 18, legendY + 70);

  // Edge legend - double line for closest
  ctx.beginPath();
  ctx.moveTo(legendX, legendY + 86);
  ctx.lineTo(legendX + 12, legendY + 86);
  ctx.moveTo(legendX, legendY + 90);
  ctx.lineTo(legendX + 12, legendY + 90);
  ctx.strokeStyle = '#999';
  ctx.lineWidth = 2;
  ctx.stroke();

  ctx.fillStyle = '#444';
  ctx.fillText('Closest to center', legendX + 18, legendY + 88);
}

// Start loading when DOM is ready
document.addEventListener('DOMContentLoaded', loadDataAndInit);

// ========== SEARCH BAR FUNCTIONALITY ==========
(function initSearchBar() {
  const searchInput = document.getElementById('pairSearch');
  const searchBtn = document.getElementById('searchBtn');
  const datalist = document.getElementById('pairSearchOptions');

  if (!searchInput || !searchBtn) return;

  // Load available pairs for autocomplete from static index
  fetch(`${DATA_BASE}/index.json`)
    .then(resp => resp.json())
    .then(pairs => {
      if (datalist && Array.isArray(pairs)) {
        // Populate datalist with first 100 pairs for performance
        pairs.slice(0, 100).forEach(p => {
          const opt = document.createElement('option');
          opt.value = p.id || p.pair_id || `${p.geneA}_${p.geneB}`;
          datalist.appendChild(opt);
        });
      }
    })
    .catch(e => console.warn('Failed to load pairs for search:', e));

  // Navigate to selected pair
  function navigateToPair() {
    const value = searchInput.value.trim();
    if (value) {
      window.location.href = `/?pair=${encodeURIComponent(value)}`;
    }
  }

  searchBtn.addEventListener('click', navigateToPair);
  searchInput.addEventListener('keydown', e => {
    if (e.key === 'Enter') navigateToPair();
  });

  // Pre-fill current pair if available
  if (PAIR_ID) {
    searchInput.value = PAIR_ID;
  }
})();

// ========== NOTEBOOK CODE BELOW ==========
// (Adapted from the original notebook - DATA, SUMMARY, PDB64_FULL are now loaded via API)


const AA_ORDER = ['K','R','H','E','D','N','Q','T','S','C','G','A','V','L','I','M','P','Y','F','W'];
let AM_MODES = ["raw"];
const MAX_SHARED_LIST = 30;
const MAX_UNIQUE_LIST = 20;
const MAX_SHARED_GRAPH = 7;
const MAX_UNIQUE_GRAPH = 4;

let amMode = 'raw';
let amTrackA = null;
let amTrackB = null;
let damTrack = null;
let amMatrixTracksA = [];
let amMatrixTracksB = [];

/* ----------------- PDBe complexes data ----------------- */
let PDBe_COMPLEXES = [];
let UNIPROT_A = "";
let UNIPROT_B = "";

/* ========== SUMMARY SECTION FUNCTIONS ========== */
let radarChart = null;
let boxplotChart = null;
let activeMetricKey = null;
let PPI_GRAPH_DATA = null;
let showUniquePpis = true;
const DEFAULT_BOXPLOT_HINT = 'Click a radar point or metric card to compare this pair with the cohort';

function initSummarySection() {
  const pair = SUMMARY.pair || {};
  const gene1 = SUMMARY.gene1 || {};
  const gene2 = SUMMARY.gene2 || {};
  
  // Update gene symbols throughout the page
  const g1 = DATA.g1 || gene1.symbol || 'Gene A';
  const g2 = DATA.g2 || gene2.symbol || 'Gene B';
  const a1 = DATA.a1 || gene1.uniprot || '';
  const a2 = DATA.a2 || gene2.uniprot || '';
  
  // Update header/title
  document.getElementById('titleMain').textContent = `${g1} vs ${g2}`;
  document.getElementById('gene1Symbol').textContent = g1;
  document.getElementById('gene2Symbol').textContent = g2;
  document.getElementById('acc1Display').textContent = a1;
  document.getElementById('acc2Display').textContent = a2;
  
  // Update legend labels
  const legendA = document.getElementById('legendA');
  const legendB = document.getElementById('legendB');
  if (legendA) legendA.textContent = g1;
  if (legendB) legendB.textContent = g2;
  
  // Update PPI labels
  const ppiLabelA = document.getElementById('ppiLabelA');
  const ppiLabelB = document.getElementById('ppiLabelB');
  if (ppiLabelA) ppiLabelA.textContent = g1;
  if (ppiLabelB) ppiLabelB.textContent = g2;
  
  // Update domain table headers
  const domAHeader = document.getElementById('domAHeader');
  const domBHeader = document.getElementById('domBHeader');
  if (domAHeader) domAHeader.textContent = `${g1} (${a1}) domain`;
  if (domBHeader) domBHeader.textContent = `${g2} (${a2}) domain`;
  
  // Update UniProt links
  const link1Up = document.getElementById('link1Up');
  const link2Up = document.getElementById('link2Up');
  const link1Pdbe = document.getElementById('link1Pdbe');
  const link2Pdbe = document.getElementById('link2Pdbe');
  if (link1Up && a1) link1Up.href = `https://www.uniprot.org/uniprotkb/${a1}`;
  if (link2Up && a2) link2Up.href = `https://www.uniprot.org/uniprotkb/${a2}`;
  if (link1Pdbe && a1) link1Pdbe.href = `https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/${a1}`;
  if (link2Pdbe && a2) link2Pdbe.href = `https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/${a2}`;
  
  if (gene1.is_essential) {
    document.getElementById('sum-essential1').textContent = 'Essential (DepMap)';
    document.getElementById('sum-essential1').className = 'essential-badge';
  }
  if (gene2.is_essential) {
    document.getElementById('sum-essential2').textContent = 'Essential (DepMap)';
    document.getElementById('sum-essential2').className = 'essential-badge';
  }
  
  if (gene1.chromosome && gene1.chromosome.chromosome && gene1.chromosome.chromosome !== 'NA') {
    const chr1 = gene1.chromosome;
    document.getElementById('chr-loc1').innerHTML = `<strong>Chr:</strong> ${chr1.chromosome} : ${Number(chr1.start).toLocaleString()} - ${Number(chr1.end).toLocaleString()}`;
  } else {
    document.getElementById('chr-loc1').textContent = 'Chromosome info not available';
  }
  
  if (gene2.chromosome && gene2.chromosome.chromosome && gene2.chromosome.chromosome !== 'NA') {
    const chr2 = gene2.chromosome;
    document.getElementById('chr-loc2').innerHTML = `<strong>Chr:</strong> ${chr2.chromosome} : ${Number(chr2.start).toLocaleString()} - ${Number(chr2.end).toLocaleString()}`;
  } else {
    document.getElementById('chr-loc2').textContent = 'Chromosome info not available';
  }
  
  const boolToStr = (v) => v === true ? 'Yes' : (v === false ? 'No' : '–');
  document.getElementById('sum-wgd').textContent = boolToStr(pair.wgd);
  document.getElementById('sum-family-size').textContent = (pair.family_size ?? '–');
  document.getElementById('sum-closest').textContent = boolToStr(pair.closest);
  document.getElementById('sum-same-chr').textContent = boolToStr(pair.same_chr);
  document.getElementById('sum-interact').textContent = boolToStr(pair.interact_bioplex);
  document.getElementById('sum-shared-ppi').textContent = (pair.n_shared_ppi ?? '0');
  
  renderPpiSection(pair, gene1, gene2);
  initRadarChart();
  renderConservationList();
  resetMetricSelection();

  // Initialize new sections
  initProteinDescriptions();
  initSimilaritySearchSection();
  initFamilyFeaturesSection();
  initSlFunctionalSection();

  const resetBtn = document.getElementById('resetMetricView');
  if (resetBtn && !resetBtn.dataset.bound) {
    resetBtn.addEventListener('click', () => resetMetricSelection(), { passive: true });
    resetBtn.dataset.bound = '1';
  }
}

function initRadarChart() {
  const conservation = SUMMARY.conservation || {};
  const wrapper = document.getElementById('radarChartWrapper');
  if (!wrapper) return false;

  const labels = [];
  const radarValues = [];
  const metricKeys = [];
  
  for (const [key, info] of Object.entries(conservation)) {
    labels.push(info.label || key);
    radarValues.push(typeof info.radar_value === "number" ? info.radar_value : 50);
    metricKeys.push(key);
  }
  
  if (labels.length === 0) {
    wrapper.innerHTML = '<div class="boxplot-hint">Conservation data not available</div>';
    return false;
  }

  if (typeof Chart === 'undefined') {
    renderStaticRadar(wrapper, labels, radarValues);
    radarChart = null;
    return false;
  }

  const canvas = ensureRadarCanvas(wrapper);
  const ctx = canvas.getContext('2d');
  if (radarChart) {
    radarChart.destroy();
    radarChart = null;
  }
  
  radarChart = new Chart(ctx, {
    type: 'radar',
    data: {
      labels: labels,
      datasets: [{
        label: 'Conservation Percentile',
        data: radarValues,
        fill: true,
        backgroundColor: 'rgba(102, 126, 234, 0.2)',
        borderColor: 'rgba(102, 126, 234, 1)',
        pointBackgroundColor: 'rgba(102, 126, 234, 1)',
        pointBorderColor: '#fff',
        pointHoverBackgroundColor: '#fff',
        pointHoverBorderColor: 'rgba(102, 126, 234, 1)',
        pointRadius: 6,
        pointHoverRadius: 8,
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      scales: {
        r: {
          angleLines: { display: true },
          suggestedMin: 0,
          suggestedMax: 100,
          ticks: { stepSize: 25, callback: (v) => v + '%' }
        }
      },
      plugins: {
        legend: { display: false },
        tooltip: { callbacks: { label: (ctx) => (((ctx.parsed && ctx.parsed.r) ?? 0).toFixed(1) + '% percentile') } }
      },
      onClick: (event, elements) => {
        if (elements.length > 0) {
          const idx = elements[0].index;
          showBoxplotForMetric(metricKeys[idx]);
        }
      }
    }
  });
  return true;
}

function ensureRadarCanvas(wrapper){
  let canvas = wrapper.querySelector('canvas#radarChart');
  if (!canvas) {
    canvas = document.createElement('canvas');
    canvas.id = 'radarChart';
    wrapper.innerHTML = '';
    wrapper.appendChild(canvas);
  }
  return canvas;
}

function renderStaticRadar(wrapper, labels, values){
  const size = 320;
  const center = size / 2;
  const radius = Math.min(center - 20, 130);
  const rings = 4;
  wrapper.innerHTML = '';

  const svgNS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(svgNS, 'svg');
  svg.setAttribute('viewBox', `0 0 ${size} ${size}`);

  const angleStep = (Math.PI * 2) / labels.length;
  const toPoint = (val, idx) => {
    const angle = -Math.PI / 2 + idx * angleStep;
    const r = (Math.max(0, Math.min(100, val)) / 100) * radius;
    return {
      x: center + r * Math.cos(angle),
      y: center + r * Math.sin(angle)
    };
  };

  // rings
  for (let i = 1; i <= rings; i++){
    const r = (i / rings) * radius;
    const circle = document.createElementNS(svgNS, 'circle');
    circle.setAttribute('cx', center);
    circle.setAttribute('cy', center);
    circle.setAttribute('r', r);
    circle.setAttribute('fill', 'none');
    circle.setAttribute('stroke', '#f0e7d3');
    circle.setAttribute('stroke-width', '1');
    svg.appendChild(circle);
  }

  // axes + labels
  labels.forEach((label, idx) => {
    const angle = -Math.PI / 2 + idx * angleStep;
    const x = center + radius * Math.cos(angle);
    const y = center + radius * Math.sin(angle);

    const axis = document.createElementNS(svgNS, 'line');
    axis.setAttribute('x1', center);
    axis.setAttribute('y1', center);
    axis.setAttribute('x2', x);
    axis.setAttribute('y2', y);
    axis.setAttribute('stroke', '#e0d5bf');
    svg.appendChild(axis);

    const text = document.createElementNS(svgNS, 'text');
    text.setAttribute('x', center + (radius + 12) * Math.cos(angle));
    text.setAttribute('y', center + (radius + 12) * Math.sin(angle));
    text.setAttribute('text-anchor', Math.cos(angle) > 0.1 ? 'start' : (Math.cos(angle) < -0.1 ? 'end' : 'middle'));
    text.setAttribute('dominant-baseline', Math.sin(angle) > 0.1 ? 'hanging' : (Math.sin(angle) < -0.1 ? 'baseline' : 'middle'));
    text.setAttribute('font-size', '11px');
    text.setAttribute('fill', '#5c4d30');
    text.textContent = label;
    svg.appendChild(text);
  });

  // polygon
  const polygon = document.createElementNS(svgNS, 'polygon');
  const pts = values.map((val, idx) => {
    const pt = toPoint(val, idx);
    return `${pt.x},${pt.y}`;
  });
  polygon.setAttribute('points', pts.join(' '));
  polygon.setAttribute('fill', 'rgba(102,126,234,0.25)');
  polygon.setAttribute('stroke', 'rgba(102,126,234,0.9)');
  polygon.setAttribute('stroke-width', '2');
  svg.appendChild(polygon);

  wrapper.appendChild(svg);
  const note = document.createElement('div');
  note.className = 'radar-fallback-note';
  note.textContent = 'Static radar preview (Chart.js unavailable).';
  wrapper.appendChild(note);
}

function renderConservationList() {
  const container = document.getElementById('conservationList');
  if (!container) return;
  const conservation = SUMMARY.conservation || {};
  const entries = Object.entries(conservation);
  if (!entries.length) {
    container.innerHTML = '<div class="boxplot-hint">Conservation metrics unavailable.</div>';
    return;
  }

  container.innerHTML = '';
  entries.forEach(([metricKey, info]) => {
    const card = document.createElement('div');
    card.className = 'metric';
    card.dataset.metric = metricKey;
    card.setAttribute('role', 'button');
    card.tabIndex = 0;
    const directionTitle = info.direction_hint || (info.higher_is_more_conserved ? 'Higher values = more conserved' : 'Lower values = more conserved');
    card.title = directionTitle;

    const value = typeof info.value === 'number' ? info.value.toFixed(3) : '–';
    const pctVal = typeof info.percentile === 'number' ? info.percentile.toFixed(1) : null;
    const pctLabel = pctVal ? `${pctVal}% percentile` : 'Percentile unavailable';
    const badgeClass = pctVal === null ? '' : (info.percentile >= 75 ? 'cons-high' : (info.percentile >= 50 ? 'cons-medium' : 'cons-low'));

    card.innerHTML = `<span class="label">${info.label || ''}</span><span class="value">${value}</span><span class="percentile ${badgeClass}">${pctLabel}</span>`;
    card.addEventListener('click', () => showBoxplotForMetric(metricKey));
    card.addEventListener('keydown', (ev) => {
      if (ev.key === 'Enter' || ev.key === ' ') {
        ev.preventDefault();
        showBoxplotForMetric(metricKey);
      }
    });
    container.appendChild(card);
  });
  highlightMetricSelection();
}

/* ========== NEW SECTION FUNCTIONS ========== */

// Charts for new sections
// (radar chart removed - replaced by PLMA alignment)

const DESC_TRUNCATE_LENGTH = 200; // Characters before truncation

function setupDescriptionToggle(funcEl, toggleBtn, fullText) {
  if (!funcEl || !toggleBtn || !fullText) return;

  if (fullText.length > DESC_TRUNCATE_LENGTH) {
    // Truncate at a word boundary
    let truncated = fullText.slice(0, DESC_TRUNCATE_LENGTH);
    const lastSpace = truncated.lastIndexOf(' ');
    if (lastSpace > DESC_TRUNCATE_LENGTH - 50) {
      truncated = truncated.slice(0, lastSpace);
    }
    truncated += '...';

    funcEl.textContent = truncated;
    funcEl.classList.add('truncated');
    funcEl.dataset.fullText = fullText;
    funcEl.dataset.truncatedText = truncated;
    toggleBtn.style.display = 'inline-block';
    toggleBtn.textContent = 'Show more';

    toggleBtn.addEventListener('click', () => {
      const isExpanded = funcEl.classList.contains('expanded');
      if (isExpanded) {
        funcEl.textContent = funcEl.dataset.truncatedText;
        funcEl.classList.remove('expanded');
        funcEl.classList.add('truncated');
        toggleBtn.textContent = 'Show more';
      } else {
        funcEl.textContent = funcEl.dataset.fullText;
        funcEl.classList.remove('truncated');
        funcEl.classList.add('expanded');
        toggleBtn.textContent = 'Show less';
      }
    }, { passive: true });
  } else {
    funcEl.textContent = fullText;
    toggleBtn.style.display = 'none';
  }
}

function initProteinDescriptions() {
  const gene1 = SUMMARY.gene1 || {};
  const gene2 = SUMMARY.gene2 || {};
  const a1 = DATA.a1 || '';
  const a2 = DATA.a2 || '';

  // Gene 1 description
  const desc1 = gene1.description || {};
  const gene1DescEl = document.getElementById('gene1Desc');
  if (gene1DescEl && (desc1.name || desc1.function)) {
    document.getElementById('gene1DescTitle').textContent = desc1.name || 'Unknown protein';

    const funcEl = document.getElementById('gene1DescFunc');
    const toggleBtn = document.getElementById('gene1DescToggle');
    const fullText = desc1.function || 'Function not available';
    setupDescriptionToggle(funcEl, toggleBtn, fullText);

    // Set UniProt source link
    const sourceLink = document.getElementById('gene1DescSource');
    if (sourceLink && a1) {
      sourceLink.href = `https://www.uniprot.org/uniprotkb/${a1}`;
    }

    gene1DescEl.style.display = 'block';
  }

  // Gene 2 description
  const desc2 = gene2.description || {};
  const gene2DescEl = document.getElementById('gene2Desc');
  if (gene2DescEl && (desc2.name || desc2.function)) {
    document.getElementById('gene2DescTitle').textContent = desc2.name || 'Unknown protein';

    const funcEl = document.getElementById('gene2DescFunc');
    const toggleBtn = document.getElementById('gene2DescToggle');
    const fullText = desc2.function || 'Function not available';
    setupDescriptionToggle(funcEl, toggleBtn, fullText);

    // Set UniProt source link
    const sourceLink = document.getElementById('gene2DescSource');
    if (sourceLink && a2) {
      sourceLink.href = `https://www.uniprot.org/uniprotkb/${a2}`;
    }

    gene2DescEl.style.display = 'block';
  }
}

// Boxplot chart instances for new sections
// (boxplot chart removed - replaced by PLMA alignment)

// Current similarity search mode (struct or seq) and view (0=overview, 1=scale)
let simSearchMode = 'struct';
let simSearchView = 0;
// Store hit-test regions for hover tooltips
let simSearchHitRegions = [];
let simSearchBarHitRegions = [];

function initSimilaritySearchSection() {
  const canvas = document.getElementById('simSearchRankCanvas');
  const barCanvas = document.getElementById('simSearchBarCanvas');
  const modeSelect = document.getElementById('simSearchModeSelect');
  if (!canvas) return;

  // Initial draw
  drawSimSearchRankViz(simSearchMode);
  drawSimSearchBarViz(simSearchMode);

  // Mode switch handler
  if (modeSelect && !modeSelect.dataset.bound) {
    modeSelect.dataset.bound = 'true';
    modeSelect.addEventListener('change', (e) => {
      simSearchMode = e.target.value;
      drawSimSearchRankViz(simSearchMode);
      drawSimSearchBarViz(simSearchMode);
    });
  }

  // View toggle buttons
  const btn0 = document.getElementById('simSearchViewBtn0');
  const btn1 = document.getElementById('simSearchViewBtn1');
  const scroll = document.getElementById('simSearchScroll');
  function setView(idx) {
    simSearchView = idx;
    if (scroll) scroll.style.transform = `translateX(-${idx * 50}%)`;
    if (btn0) btn0.style.background = idx === 0 ? '#e8e4ff' : '';
    if (btn0) btn0.style.color = idx === 0 ? '#4f46e5' : '';
    if (btn0) btn0.style.fontWeight = idx === 0 ? '600' : '';
    if (btn1) btn1.style.background = idx === 1 ? '#e8e4ff' : '';
    if (btn1) btn1.style.color = idx === 1 ? '#4f46e5' : '';
    if (btn1) btn1.style.fontWeight = idx === 1 ? '600' : '';
  }
  if (btn0) btn0.addEventListener('click', () => setView(0));
  if (btn1) btn1.addEventListener('click', () => setView(1));
  setView(0);

  // Hover tooltip for overview canvas
  const tooltip = document.getElementById('simSearchTooltip');
  if (canvas && tooltip) {
    canvas.addEventListener('mousemove', (e) => {
      const rect = canvas.getBoundingClientRect();
      const x = (e.clientX - rect.left) * (canvas.width / (rect.width * (Math.max(2, window.devicePixelRatio || 1))));
      const y = (e.clientY - rect.top) * (canvas.height / (rect.height * (Math.max(2, window.devicePixelRatio || 1))));
      let hit = null;
      for (const region of simSearchHitRegions) {
        if (x >= region.x && x <= region.x + region.w && y >= region.y && y <= region.y + region.h) {
          hit = region;
          break;
        }
      }
      if (hit) {
        tooltip.innerHTML = hit.tooltip;
        tooltip.style.display = 'block';
        tooltip.style.left = (e.clientX + 12) + 'px';
        tooltip.style.top = (e.clientY - 8) + 'px';
        tooltip.style.position = 'fixed';
        canvas.style.cursor = 'default';
      } else {
        tooltip.style.display = 'none';
        canvas.style.cursor = 'default';
      }
    });
    canvas.addEventListener('mouseleave', () => { tooltip.style.display = 'none'; });
  }

  // Hover tooltip for scale bar canvas
  if (barCanvas && tooltip) {
    barCanvas.addEventListener('mousemove', (e) => {
      const rect = barCanvas.getBoundingClientRect();
      const dpr = Math.max(2, window.devicePixelRatio || 1);
      const x = (e.clientX - rect.left) * (barCanvas.width / (rect.width * dpr));
      const y = (e.clientY - rect.top) * (barCanvas.height / (rect.height * dpr));
      let hit = null;
      for (const region of simSearchBarHitRegions) {
        if (x >= region.x && x <= region.x + region.w && y >= region.y && y <= region.y + region.h) {
          hit = region;
          break;
        }
      }
      if (hit) {
        tooltip.innerHTML = hit.tooltip;
        tooltip.style.display = 'block';
        tooltip.style.left = (e.clientX + 12) + 'px';
        tooltip.style.top = (e.clientY - 8) + 'px';
        tooltip.style.position = 'fixed';
        barCanvas.style.cursor = 'default';
      } else {
        tooltip.style.display = 'none';
        barCanvas.style.cursor = 'default';
      }
    });
    barCanvas.addEventListener('mouseleave', () => { tooltip.style.display = 'none'; });
  }
}

// Helper: format number with commas
function formatNum(n) {
  if (n == null) return '-';
  return Number(n).toLocaleString('en-US');
}

// Shared helpers for both canvases
function ssRoundRect(ctx, x, y, w, h, r) {
  if (w <= 0) return;
  r = Math.min(r, w / 2, h / 2);
  ctx.beginPath();
  ctx.moveTo(x + r, y);
  ctx.lineTo(x + w - r, y);
  ctx.quadraticCurveTo(x + w, y, x + w, y + r);
  ctx.lineTo(x + w, y + h - r);
  ctx.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
  ctx.lineTo(x + r, y + h);
  ctx.quadraticCurveTo(x, y + h, x, y + h - r);
  ctx.lineTo(x, y + r);
  ctx.quadraticCurveTo(x, y, x + r, y);
  ctx.closePath();
}

function ssGetMetrics(mode) {
  const simSearch = SUMMARY?.similarity_search || {};
  const suffix = mode === 'struct' ? '_struct' : '_seq';
  return {
    suffix,
    dbFullName: mode === 'struct' ? 'AlphaFold DB' : 'UniProt',
    rankInfo: simSearch['rank' + suffix],
    selfSPInfo: simSearch['selfSP' + suffix],
    taxidInfo: simSearch['taxid' + suffix],
    gene1: SUMMARY?.gene1?.symbol || 'A',
    gene2: SUMMARY?.gene2?.symbol || 'B',
  };
}

// Compute bar fill % for rank: linear scale capped at p95
function rankLinearPct(value, p95) {
  if (value == null || p95 == null || p95 <= 0) return 0;
  return Math.min(100, (value / p95) * 100);
}

// ====== VIEW 1: Overview (A - DB - B) ======
function drawSimSearchRankViz(mode) {
  const canvas = document.getElementById('simSearchRankCanvas');
  const legendEl = document.getElementById('simSearchLegend');
  if (!canvas || !SUMMARY) return;

  const ctx = canvas.getContext('2d');
  const dpr = Math.max(2, window.devicePixelRatio || 1);
  simSearchHitRegions = [];

  const displayWidth = 680;
  const displayHeight = 280;
  canvas.width = displayWidth * dpr;
  canvas.height = displayHeight * dpr;
  canvas.style.width = displayWidth + 'px';
  canvas.style.height = displayHeight + 'px';
  ctx.scale(dpr, dpr);

  // Background matching page
  ctx.fillStyle = '#fafafa';
  ctx.fillRect(0, 0, displayWidth, displayHeight);

  const m = ssGetMetrics(mode);
  const rank = m.rankInfo?.value ?? null;
  const selfSP = m.selfSPInfo?.value ?? null;
  const taxid = m.taxidInfo?.value ?? null;

  // Linear scale for rank (value / p95_cap)
  const rankP95 = m.rankInfo?.p95_value ?? 500;
  const rankFillPct = rankLinearPct(rank, rankP95);

  // Percentiles for selfSP and taxid
  const selfSPPct = m.selfSPInfo?.percentile ?? 50;
  const selfSPPctRel = m.selfSPInfo?.percentile_rank_relative ?? selfSPPct;
  const taxidPct = m.taxidInfo?.percentile ?? 50;
  const taxidPctRel = m.taxidInfo?.percentile_rank_relative ?? taxidPct;

  // Total pairs for tooltip
  const totalPairs = m.rankInfo?.total_pairs ?? 105107;

  // PPI network-unified colors
  const geneAColor = '#d97706';   // Amber (PPI center gene)
  const geneAStroke = '#92400e';
  const geneBColor = '#7c3aed';   // Purple (PPI partner)
  const geneBStroke = '#5b21b6';
  // Page-matching box colors (cream/grey, no coral)
  const dbBoxColor = '#f5f1e6';
  const dbBoxStroke = '#e0d7c2';
  const barBgColor = '#eae6dc';
  const barFillColor = '#5f4d2f';
  const barFillRelColor = '#8b7a5e';
  const barBorderColor = '#333';

  // Layout
  const centerY = displayHeight / 2;
  const geneAx = 80;
  const geneBx = displayWidth - 80;
  const dbCenterX = displayWidth / 2;
  const geneRadius = 14;  // PPI network size

  // Database box
  const dbBoxWidth = 260;
  const dbBoxHeight = 210;
  const dbBoxLeft = dbCenterX - dbBoxWidth / 2;
  const dbBoxTop = centerY - dbBoxHeight / 2;
  const dbBoxRadius = 12;

  // Helper: draw bar with black border
  function drawBar(x, y, width, height, fillPct, fillColor, r) {
    r = r || 5;
    // Full background
    ssRoundRect(ctx, x, y, width, height, r);
    ctx.fillStyle = barBgColor;
    ctx.fill();
    ctx.strokeStyle = barBorderColor;
    ctx.lineWidth = 1;
    ctx.stroke();
    // Fill
    if (fillPct > 0) {
      const fw = Math.max(2, (fillPct / 100) * width);
      ssRoundRect(ctx, x, y, fw, height, r);
      ctx.fillStyle = fillColor;
      ctx.fill();
      ctx.strokeStyle = barBorderColor;
      ctx.lineWidth = 1;
      ctx.stroke();
    }
  }

  // Helper: draw split bar (top + bottom with gap)
  function drawSplitBar(x, y, width, pctTop, pctBottom) {
    const bh = 12;
    const gap = 3;
    drawBar(x, y, width, bh, pctTop, barFillColor, 4);
    drawBar(x, y + bh + gap, width, bh, pctBottom, barFillRelColor, 4);
  }

  // Gene A circle (PPI style: glow + circle + label below)
  ctx.beginPath();
  ctx.arc(geneAx, centerY, geneRadius + 6, 0, Math.PI * 2);
  ctx.fillStyle = geneAColor + '30';
  ctx.fill();
  ctx.beginPath();
  ctx.arc(geneAx, centerY, geneRadius, 0, Math.PI * 2);
  ctx.fillStyle = geneAColor;
  ctx.fill();
  ctx.strokeStyle = geneAStroke;
  ctx.lineWidth = 2;
  ctx.stroke();
  // Label below circle
  ctx.font = 'bold 12px -apple-system, sans-serif';
  ctx.fillStyle = geneAStroke;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText(m.gene1, geneAx, centerY + geneRadius + 8);

  // Gene B circle (PPI style)
  ctx.beginPath();
  ctx.arc(geneBx, centerY, geneRadius + 6, 0, Math.PI * 2);
  ctx.fillStyle = geneBColor + '30';
  ctx.fill();
  ctx.beginPath();
  ctx.arc(geneBx, centerY, geneRadius, 0, Math.PI * 2);
  ctx.fillStyle = geneBColor;
  ctx.fill();
  ctx.strokeStyle = geneBStroke;
  ctx.lineWidth = 2;
  ctx.stroke();
  ctx.font = 'bold 12px -apple-system, sans-serif';
  ctx.fillStyle = geneBStroke;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText(m.gene2, geneBx, centerY + geneRadius + 8);

  // Connecting lines A <--- DB ---> B
  ctx.strokeStyle = '#bbb';
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(geneAx + geneRadius + 8, centerY);
  ctx.lineTo(dbBoxLeft - 6, centerY);
  ctx.stroke();
  ctx.beginPath();
  ctx.moveTo(dbBoxLeft + dbBoxWidth + 6, centerY);
  ctx.lineTo(geneBx - geneRadius - 8, centerY);
  ctx.stroke();

  // Database box (cream, no coral)
  ssRoundRect(ctx, dbBoxLeft, dbBoxTop, dbBoxWidth, dbBoxHeight, dbBoxRadius);
  ctx.fillStyle = dbBoxColor;
  ctx.fill();
  ctx.strokeStyle = dbBoxStroke;
  ctx.lineWidth = 1.5;
  ctx.stroke();

  // DB name at top
  ctx.font = 'bold 14px -apple-system, sans-serif';
  ctx.fillStyle = '#5f4d2f';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText(m.dbFullName, dbCenterX, dbBoxTop + 12);

  // Metrics inside DB box
  const metricStartY = dbBoxTop + 42;
  const metricSpacing = 56;
  const labelX = dbBoxLeft + 14;
  const barX = dbBoxLeft + 14;
  const barWidth = dbBoxWidth - 28;

  const metricLabels = [
    'Rank of the paralog',
    'Human proteins ranking better',
    'Species with proteins ranking better',
  ];
  const metricValues = [rank, selfSP, taxid];

  // --- Rank row ---
  const ry = metricStartY;
  ctx.font = '12px -apple-system, sans-serif';
  ctx.fillStyle = '#666';
  ctx.textAlign = 'left';
  ctx.textBaseline = 'middle';
  ctx.fillText(metricLabels[0], labelX, ry);
  ctx.font = 'bold 14px -apple-system, sans-serif';
  ctx.fillStyle = '#333';
  ctx.textAlign = 'right';
  ctx.fillText(rank !== null ? formatNum(rank) : '-', barX + barWidth, ry);
  if (rank !== null) {
    drawBar(barX, ry + 8, barWidth, 14, rankFillPct, barFillColor, 5);
    const rPos = m.rankInfo?.rank_position ?? null;
    simSearchHitRegions.push({
      x: barX, y: ry + 8, w: barWidth, h: 14,
      tooltip: rPos != null ? `Rank ${formatNum(rank)}: ${formatNum(rPos)} of ${formatNum(totalPairs)} paralog pairs have a rank ≤ ${formatNum(rank)}` : `Rank: ${formatNum(rank)}`,
    });
  }

  // --- selfSP row (split bar) ---
  const sy = metricStartY + metricSpacing;
  ctx.font = '12px -apple-system, sans-serif';
  ctx.fillStyle = '#666';
  ctx.textAlign = 'left';
  ctx.fillText(metricLabels[1], labelX, sy);
  ctx.font = 'bold 14px -apple-system, sans-serif';
  ctx.fillStyle = '#333';
  ctx.textAlign = 'right';
  ctx.fillText(selfSP !== null ? formatNum(selfSP) : '-', barX + barWidth, sy);
  if (selfSP !== null) {
    drawSplitBar(barX, sy + 8, barWidth, selfSPPct, selfSPPctRel);
    const sPos = m.selfSPInfo?.rank_position ?? null;
    simSearchHitRegions.push({
      x: barX, y: sy + 8, w: barWidth, h: 27,
      tooltip: sPos != null
        ? `${formatNum(selfSP)} human proteins rank better: ${formatNum(sPos)} of ${formatNum(totalPairs)} pairs have ≤ ${formatNum(selfSP)}<br><small>Top: vs all pairs · Bottom: vs pairs with similar rank</small>`
        : `selfSP: ${formatNum(selfSP)}`,
    });
  }

  // --- taxid row (split bar) ---
  const ty = metricStartY + metricSpacing * 2;
  ctx.font = '12px -apple-system, sans-serif';
  ctx.fillStyle = '#666';
  ctx.textAlign = 'left';
  ctx.fillText(metricLabels[2], labelX, ty);
  ctx.font = 'bold 14px -apple-system, sans-serif';
  ctx.fillStyle = '#333';
  ctx.textAlign = 'right';
  ctx.fillText(taxid !== null ? formatNum(taxid) : '-', barX + barWidth, ty);
  if (taxid !== null) {
    drawSplitBar(barX, ty + 8, barWidth, taxidPct, taxidPctRel);
    const tPos = m.taxidInfo?.rank_position ?? null;
    simSearchHitRegions.push({
      x: barX, y: ty + 8, w: barWidth, h: 27,
      tooltip: tPos != null
        ? `${formatNum(taxid)} species rank better: ${formatNum(tPos)} of ${formatNum(totalPairs)} pairs have ≤ ${formatNum(taxid)}<br><small>Top: vs all pairs · Bottom: vs pairs with similar rank</small>`
        : `taxid: ${formatNum(taxid)}`,
    });
  }

  // Update legend
  if (legendEl) {
    let legendHTML = '';
    if (rank !== null) {
      if (rank === 0) {
        legendHTML = `<strong style="color:#2e7d32">Direct match!</strong> ${m.gene2} is the top hit when searching ${m.dbFullName} with ${m.gene1}.`;
      } else {
        legendHTML = `<strong>${formatNum(rank)}</strong> proteins between ${m.gene1} and ${m.gene2}.`;
        if (selfSP !== null) legendHTML += ` <strong>${formatNum(selfSP)}</strong> human.`;
        if (taxid !== null) legendHTML += ` <strong>${formatNum(taxid)}</strong> species.`;
        legendHTML += '<br><span style="color:#888">Bars: empty = close paralogs, filled = distant. Split bars: top = vs all pairs, bottom = vs similar rank.</span>';
      }
    } else {
      legendHTML = 'Rank data not available for this pair.';
    }
    legendEl.innerHTML = legendHTML;
  }
}

// ====== VIEW 2: Scale bars (A====B - - - - z) ======
function drawSimSearchBarViz(mode) {
  const canvas = document.getElementById('simSearchBarCanvas');
  if (!canvas || !SUMMARY) return;

  const ctx = canvas.getContext('2d');
  const dpr = Math.max(2, window.devicePixelRatio || 1);
  simSearchBarHitRegions = [];

  const displayWidth = 680;
  const displayHeight = 220;
  canvas.width = displayWidth * dpr;
  canvas.height = displayHeight * dpr;
  canvas.style.width = displayWidth + 'px';
  canvas.style.height = displayHeight + 'px';
  ctx.scale(dpr, dpr);

  ctx.fillStyle = '#fafafa';
  ctx.fillRect(0, 0, displayWidth, displayHeight);

  const m = ssGetMetrics(mode);
  const rank = m.rankInfo?.value ?? null;
  const selfSP = m.selfSPInfo?.value ?? null;
  const taxid = m.taxidInfo?.value ?? null;

  // Use max values from data for the full bar width
  const rankMax = m.rankInfo?.max_value ?? 1;
  const selfSPMax = m.selfSPInfo?.max_value ?? 1;
  const taxidMax = m.taxidInfo?.max_value ?? 1;
  const totalPairs = m.rankInfo?.total_pairs ?? 105107;

  // PPI colors
  const geneAColor = '#d97706';
  const geneAStroke = '#92400e';
  const geneBColor = '#7c3aed';
  const geneBStroke = '#5b21b6';
  const filledColor = '#e8dcc8';
  const filledBorder = '#333';
  const emptyBorder = '#999';

  const barLabels = [
    'Rank of the paralog',
    'Human proteins ranking better',
    'Species with proteins ranking better',
  ];
  const values = [rank, selfSP, taxid];
  const maxValues = [rankMax, selfSPMax, taxidMax];
  const infos = [m.rankInfo, m.selfSPInfo, m.taxidInfo];

  const padLeft = 24;
  const padRight = 24;
  const barTotalWidth = displayWidth - padLeft - padRight;
  const barHeight = 28;
  const barSpacing = 62;
  const startY = 28;
  const circleR = 12;

  for (let i = 0; i < 3; i++) {
    const val = values[i];
    const maxVal = maxValues[i];
    const info = infos[i];
    const by = startY + i * barSpacing;

    // Label above bar
    ctx.font = '12px -apple-system, sans-serif';
    ctx.fillStyle = '#666';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'bottom';
    ctx.fillText(barLabels[i], padLeft, by - 3);
    // Value text right-aligned
    ctx.textAlign = 'right';
    ctx.fillStyle = '#888';
    ctx.fillText(val != null ? `${formatNum(val)} / ${formatNum(maxVal)}` : '-', displayWidth - padRight, by - 3);

    if (val == null || maxVal <= 0) continue;

    const fillFrac = Math.min(1, val / maxVal);
    const fillW = Math.max(0, fillFrac * barTotalWidth);
    const emptyW = barTotalWidth - fillW;

    // Filled portion (A====B)
    if (fillW > 0) {
      ssRoundRect(ctx, padLeft, by, fillW, barHeight, 5);
      ctx.fillStyle = filledColor;
      ctx.fill();
      ctx.strokeStyle = filledBorder;
      ctx.lineWidth = 1.5;
      ctx.stroke();
    }

    // Empty portion (- - - z) with dashed border
    if (emptyW > 2) {
      ssRoundRect(ctx, padLeft + fillW, by, emptyW, barHeight, 5);
      ctx.fillStyle = '#fafafa';
      ctx.fill();
      ctx.setLineDash([4, 4]);
      ctx.strokeStyle = emptyBorder;
      ctx.lineWidth = 1;
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // Gene A circle at start
    const ax = padLeft;
    const cy = by + barHeight / 2;
    ctx.beginPath();
    ctx.arc(ax, cy, circleR, 0, Math.PI * 2);
    ctx.fillStyle = geneAColor;
    ctx.fill();
    ctx.strokeStyle = geneAStroke;
    ctx.lineWidth = 1.5;
    ctx.stroke();
    ctx.font = 'bold 9px -apple-system, sans-serif';
    ctx.fillStyle = '#fff';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(m.gene1, ax, cy);

    // Gene B circle at fill boundary
    const bx = padLeft + fillW;
    ctx.beginPath();
    ctx.arc(bx, cy, circleR, 0, Math.PI * 2);
    ctx.fillStyle = geneBColor;
    ctx.fill();
    ctx.strokeStyle = geneBStroke;
    ctx.lineWidth = 1.5;
    ctx.stroke();
    ctx.font = 'bold 9px -apple-system, sans-serif';
    ctx.fillStyle = '#fff';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(m.gene2, bx, cy);

    // "max" label at end
    ctx.font = '10px -apple-system, sans-serif';
    ctx.fillStyle = '#999';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    ctx.fillText('max', displayWidth - padRight - 4, cy);

    // Hit region for tooltip
    const rPos = info?.rank_position ?? null;
    const pctStr = (fillFrac * 100).toFixed(1);
    simSearchBarHitRegions.push({
      x: padLeft, y: by, w: barTotalWidth, h: barHeight,
      tooltip: rPos != null
        ? `${formatNum(val)} out of max ${formatNum(maxVal)} (${pctStr}%)<br>${formatNum(rPos)} of ${formatNum(totalPairs)} pairs have ≤ ${formatNum(val)}`
        : `${formatNum(val)} out of max ${formatNum(maxVal)} (${pctStr}%)`,
    });
  }
}

// ========== PLMA ALIGNMENT VISUALIZATION ==========
let plmaHitRegions = [];

function initFamilyFeaturesSection() {
  const canvas = document.getElementById('plmaAlignCanvas');
  if (!canvas || !PLMA_DATA) {
    const body = document.getElementById('familyFeaturesBody');
    if (body && !PLMA_DATA) {
      body.innerHTML = '<p class="small" style="color:#888;">PLMA alignment data not available for this pair.</p>';
    }
    return;
  }

  drawPlmaAlignment();

  // Tooltip handling
  const tooltip = document.getElementById('plmaTooltip');
  if (tooltip) {
    canvas.addEventListener('mousemove', (e) => {
      const rect = canvas.getBoundingClientRect();
      const dpr = Math.max(2, window.devicePixelRatio || 1);
      const x = (e.clientX - rect.left) * (canvas.width / (rect.width * dpr));
      const y = (e.clientY - rect.top) * (canvas.height / (rect.height * dpr));
      let hit = null;
      for (const region of plmaHitRegions) {
        if (x >= region.x && x <= region.x + region.w && y >= region.y && y <= region.y + region.h) {
          hit = region;
          break;
        }
      }
      if (hit) {
        tooltip.innerHTML = hit.tooltip;
        tooltip.style.display = 'block';
        tooltip.style.left = (e.clientX + 14) + 'px';
        tooltip.style.top = (e.clientY - 6) + 'px';
        canvas.style.cursor = 'default';
      } else {
        tooltip.style.display = 'none';
      }
    });
    canvas.addEventListener('mouseleave', () => { tooltip.style.display = 'none'; });
  }
}

function drawPlmaAlignment() {
  const canvas = document.getElementById('plmaAlignCanvas');
  const legendEl = document.getElementById('plmaLegend');
  const summaryEl = document.getElementById('plmaSummary');
  if (!canvas || !PLMA_DATA) return;

  const ctx = canvas.getContext('2d');
  const dpr = Math.max(2, window.devicePixelRatio || 1);
  plmaHitRegions = [];

  const plma = PLMA_DATA;
  const sequences = plma.sequences || [];
  const blocks = plma.blocks || [];
  const geneASeq = plma.gene_a_seq;
  const geneBSeq = plma.gene_b_seq;
  const geneA = plma.gene_a;
  const geneB = plma.gene_b;

  // Reorder: gene A first, gene B second, then others
  const orderedSeqs = [];
  const seqA = sequences.find(s => s.num === geneASeq);
  const seqB = sequences.find(s => s.num === geneBSeq);
  if (seqA) orderedSeqs.push(seqA);
  if (seqB) orderedSeqs.push(seqB);
  for (const s of sequences) {
    if (s.num !== geneASeq && s.num !== geneBSeq) orderedSeqs.push(s);
  }

  const nSeqs = orderedSeqs.length;
  const nBlocks = blocks.length;

  // === Column-based layout ===
  // Each block gets a column proportional to its max length.
  // Between blocks: a gap column with dashed connectors.
  const labelWidth = 90;
  const padRight = 16;
  const padTop = 20;
  const trackHeight = nSeqs <= 6 ? 22 : (nSeqs <= 15 ? 14 : 10);
  const trackGap = nSeqs <= 6 ? 10 : (nSeqs <= 15 ? 6 : 3);
  const pairTrackHeight = trackHeight + 6;

  const scrollWrapper = document.getElementById('plmaScrollWrapper');
  const containerWidth = Math.max(700, scrollWrapper?.clientWidth || canvas.parentElement?.clientWidth || 700);

  // Compute max length per block (across all seqs in that block)
  const blockMaxLen = blocks.map(b => {
    let mx = 0;
    for (const pos of Object.values(b.positions)) mx = Math.max(mx, pos.length || 0);
    return mx;
  });
  const totalBlockAA = blockMaxLen.reduce((a, b) => a + b, 0) || 1;

  // Ensure a minimum pixel-per-AA so blocks remain readable.
  // If needed, the canvas expands beyond the container (horizontal scroll).
  const nGaps = Math.max(0, nBlocks - 1);
  const gapWidthBase = nBlocks <= 10 ? 14 : (nBlocks <= 50 ? 8 : 4);
  const minPxPerAA = 0.6;
  const minBlockArea = totalBlockAA * minPxPerAA;
  const minNeededWidth = labelWidth + padRight + minBlockArea + nGaps * gapWidthBase;
  const displayWidth = Math.max(containerWidth, minNeededWidth);
  const trackAreaWidth = displayWidth - labelWidth - padRight;

  const totalGapWidth = nGaps * gapWidthBase;
  const blockAreaWidth = Math.max(trackAreaWidth - totalGapWidth, trackAreaWidth * 0.5);
  const gapWidth = nGaps > 0 ? (trackAreaWidth - blockAreaWidth) / nGaps : 0;
  const minBlockPx = 3;

  // Column x positions (start of each block column)
  const blockColX = [];
  const blockColW = [];
  let curX = labelWidth;
  for (let bi = 0; bi < nBlocks; bi++) {
    blockColX.push(curX);
    const w = Math.max(minBlockPx, (blockMaxLen[bi] / totalBlockAA) * blockAreaWidth);
    blockColW.push(w);
    curX += w;
    if (bi < nBlocks - 1) curX += gapWidth; // gap between blocks
  }

  // Canvas height
  const totalTrackHeight = pairTrackHeight * 2 + (trackGap * 2) +
    Math.max(0, nSeqs - 2) * (trackHeight + trackGap) + padTop + 30;
  const displayHeight = Math.max(180, totalTrackHeight);

  canvas.width = displayWidth * dpr;
  canvas.height = displayHeight * dpr;
  canvas.style.width = displayWidth + 'px';
  canvas.style.height = displayHeight + 'px';
  // Ensure the inner wrapper expands so the scroll container works
  const inner = document.getElementById('plmaCanvasInner');
  if (inner) inner.style.minWidth = displayWidth + 'px';
  ctx.scale(dpr, dpr);

  ctx.fillStyle = '#fafafa';
  ctx.fillRect(0, 0, displayWidth, displayHeight);

  // Category colors
  const catColors = {
    shared_with_family: '#8b6f3a',
    pair_exclusive:     '#43a047',
    specific_a:         '#d97706',
    specific_b:         '#7c3aed',
    family_only:        '#d4cfc5',
  };
  const catBorders = {
    shared_with_family: '#5f4d2f',
    pair_exclusive:     '#2e7d32',
    specific_a:         '#92400e',
    specific_b:         '#5b21b6',
    family_only:        '#a8a298',
  };
  const catLabels = {
    shared_with_family: 'Shared (pair + family)',
    pair_exclusive:     'Pair exclusive',
    specific_a:         `Specific to ${geneA}`,
    specific_b:         `Specific to ${geneB}`,
    family_only:        'Other family members only',
  };

  // For each sequence, build a list of which block indices it participates in
  const seqBlockIndices = orderedSeqs.map(seq =>
    blocks.map((b, i) => b.positions[seq.num] ? i : -1).filter(i => i >= 0)
  );

  // Draw tracks
  let yPos = padTop;

  for (let si = 0; si < orderedSeqs.length; si++) {
    const seq = orderedSeqs[si];
    const isPairA = seq.num === geneASeq;
    const isPairB = seq.num === geneBSeq;
    const isPair = isPairA || isPairB;
    const th = isPair ? pairTrackHeight : trackHeight;
    const myBlockIndices = seqBlockIndices[si];

    // Track label
    ctx.font = isPair ? 'bold 12px -apple-system, sans-serif' : '11px -apple-system, sans-serif';
    ctx.fillStyle = isPairA ? '#92400e' : (isPairB ? '#5b21b6' : '#666');
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    const label = seq.gene || seq.uniprot || `Seq ${seq.num}`;
    ctx.fillText(label, labelWidth - 8, yPos + th / 2);

    // Draw gap connectors (dashed lines between consecutive blocks this seq has)
    const cy = yPos + th / 2;
    ctx.setLineDash([3, 3]);
    ctx.strokeStyle = isPair ? '#c0b69e' : '#d5d5d5';
    ctx.lineWidth = isPair ? 1.2 : 0.8;

    for (let k = 0; k < myBlockIndices.length - 1; k++) {
      const biLeft = myBlockIndices[k];
      const biRight = myBlockIndices[k + 1];
      const x1 = blockColX[biLeft] + blockColW[biLeft]; // right edge of left block
      const x2 = blockColX[biRight]; // left edge of right block
      if (x2 > x1 + 1) {
        ctx.beginPath();
        ctx.moveTo(x1, cy);
        ctx.lineTo(x2, cy);
        ctx.stroke();
      }
    }
    ctx.setLineDash([]);

    // Draw blocks at their aligned column positions
    for (const bi of myBlockIndices) {
      const block = blocks[bi];
      const pos = block.positions[seq.num];
      const cat = block.category;
      const bx = blockColX[bi];
      const bw = blockColW[bi];

      ctx.fillStyle = catColors[cat] || '#ccc';
      plmaRoundRect(ctx, bx, yPos + 1, bw, th - 2, 2);
      ctx.fill();
      ctx.strokeStyle = catBorders[cat] || '#999';
      ctx.lineWidth = isPair ? 1.2 : 0.8;
      ctx.stroke();

      // AA sequence for tooltip (wrap long sequences)
      const aaSeq = pos.seq || '';
      let aaHtml = '';
      if (aaSeq.length > 0) {
        // Wrap every 40 chars
        const wrapped = aaSeq.match(/.{1,40}/g) || [aaSeq];
        aaHtml = `<code style="font-size:10px;color:#555;word-break:break-all;line-height:1.3;display:block;margin-top:3px;">${wrapped.join('<br>')}</code>`;
      }

      plmaHitRegions.push({
        x: bx, y: yPos, w: bw, h: th,
        tooltip: `<strong>${block.id}</strong> · ${catLabels[cat] || cat}<br>`
          + `${label}: pos ${pos.start}–${pos.end} (${pos.length} aa)<br>`
          + `<span style="color:#888">${block.n_seqs} of ${nSeqs} family members</span>`
          + aaHtml,
      });
    }

    // Separator after pair tracks
    if (si === 1 && nSeqs > 2) {
      yPos += th + trackGap;
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(labelWidth, yPos - trackGap / 2);
      ctx.lineTo(displayWidth - padRight, yPos - trackGap / 2);
      ctx.strokeStyle = '#ddd';
      ctx.lineWidth = 1;
      ctx.stroke();
      ctx.setLineDash([]);
    } else {
      yPos += th + trackGap;
    }
  }

  // Legend
  if (legendEl) {
    const usedCats = new Set(blocks.map(b => b.category));
    let html = '';
    for (const [cat, label] of Object.entries(catLabels)) {
      if (!usedCats.has(cat)) continue;
      html += `<span style="display:inline-flex;align-items:center;gap:4px;">`
        + `<span style="display:inline-block;width:14px;height:10px;border-radius:2px;background:${catColors[cat]};border:1px solid ${catBorders[cat]}"></span>`
        + `<span>${label}</span></span>`;
    }
    // Add gap connector to legend
    html += `<span style="display:inline-flex;align-items:center;gap:4px;">`
      + `<span style="display:inline-block;width:14px;border-top:1.5px dashed #b0a890;height:0;"></span>`
      + `<span>Gap between blocks</span></span>`;
    legendEl.innerHTML = html;
  }

  // Summary
  if (summaryEl) {
    const s = plma.summary || {};
    const parts = [];
    if (s.shared_with_family) parts.push(`${s.shared_with_family} aa shared with family`);
    if (s.pair_exclusive) parts.push(`${s.pair_exclusive} aa pair-exclusive`);
    if (s.specific_a) parts.push(`${s.specific_a} aa specific to ${geneA}`);
    if (s.specific_b) parts.push(`${s.specific_b} aa specific to ${geneB}`);
    if (s.family_only) parts.push(`${s.family_only} aa in other family members only`);
    summaryEl.textContent = `${blocks.length} conserved blocks across ${nSeqs} family members` +
      (parts.length ? ` · ${parts.join(' · ')}` : '');
  }
}

function plmaRoundRect(ctx, x, y, w, h, r) {
  if (w <= 0) return;
  r = Math.min(r, w / 2, h / 2);
  ctx.beginPath();
  ctx.moveTo(x + r, y);
  ctx.lineTo(x + w - r, y);
  ctx.quadraticCurveTo(x + w, y, x + w, y + r);
  ctx.lineTo(x + w, y + h - r);
  ctx.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
  ctx.lineTo(x + r, y + h);
  ctx.quadraticCurveTo(x, y + h, x, y + h - r);
  ctx.lineTo(x, y + r);
  ctx.quadraticCurveTo(x, y, x + r, y);
  ctx.closePath();
}

// Generic boxplot drawing function for reuse
function drawGenericBoxplot(containerId, metricInfo, boxplotData, existingChart, setChart) {
  const container = document.getElementById(containerId);
  container.innerHTML = '<canvas style="width:100%;height:180px;"></canvas>';

  const canvas = container.querySelector('canvas');
  if (!canvas || typeof Chart === 'undefined') {
    container.innerHTML = '<div class="boxplot-hint">Chart.js required</div>';
    return;
  }

  const ctx = canvas.getContext('2d');
  const {q1, median, q3, whisker_low, whisker_high, pair_value} = boxplotData;

  if (existingChart) existingChart.destroy();

  const chart = new Chart(ctx, {
    type: 'bar',
    data: {
      labels: ['Distribution'],
      datasets: [
        { label: 'Lower', data: [q1 - whisker_low], backgroundColor: 'rgba(200,200,200,0.3)', barPercentage: 0.5 },
        { label: 'Q1-Med', data: [median - q1], backgroundColor: 'rgba(102, 126, 234, 0.4)', barPercentage: 0.5 },
        { label: 'Med-Q3', data: [q3 - median], backgroundColor: 'rgba(118, 75, 162, 0.4)', barPercentage: 0.5 },
        { label: 'Upper', data: [whisker_high - q3], backgroundColor: 'rgba(200,200,200,0.3)', barPercentage: 0.5 },
      ]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      indexAxis: 'y',
      scales: {
        x: { stacked: true, min: whisker_low - (whisker_high - whisker_low) * 0.1, max: whisker_high + (whisker_high - whisker_low) * 0.1, title: { display: true, text: metricInfo.label } },
        y: { stacked: true, display: false }
      },
      plugins: { legend: { display: false }, tooltip: { enabled: false } }
    },
    plugins: [{
      id: 'pairMarker',
      afterDraw: (chart) => {
        if (pair_value == null) return;
        const ctx = chart.ctx;
        const xAxis = chart.scales.x;
        const yAxis = chart.scales.y;
        const x = xAxis.getPixelForValue(pair_value);
        const y = yAxis.getPixelForValue(0);
        ctx.save();
        ctx.beginPath();
        ctx.arc(x, y, 10, 0, Math.PI * 2);
        ctx.fillStyle = '#ef5350';
        ctx.fill();
        ctx.strokeStyle = '#c62828';
        ctx.lineWidth = 2;
        ctx.stroke();
        ctx.restore();
      }
    }]
  });

  setChart(chart);
}

function initSlFunctionalSection() {
  const slFunc = SUMMARY.sl_functional || {};
  const isSL = slFunc.is_sl;
  const slScreens = slFunc.sl_screens || [];
  const goSim = slFunc.go_similarity || {};

  // SL Status - single flag with screens
  const flagsContainer = document.getElementById('slFlags');
  if (flagsContainer) {
    let html = '';

    if (isSL === true) {
      html = '<span class="sl-flag positive">Synthetic Lethal</span>';
      if (slScreens.length > 0) {
        html += `<span class="sl-note">Found in: ${slScreens.join(', ')}</span>`;
      }
    } else if (isSL === false) {
      html = '<span class="sl-flag negative">Not Synthetic Lethal</span>';
    } else {
      html = '<span class="sl-flag unknown">SL status unknown</span>';
    }

    flagsContainer.innerHTML = html;
  }

  // Hide screens container (now shown inline with flag)
  const screensContainer = document.getElementById('slScreens');
  if (screensContainer) {
    screensContainer.style.display = 'none';
  }

  // GO Similarity with color coding
  const formatGO = (key, valueElId, pctElId) => {
    const info = goSim[key] || {};
    const valueEl = document.getElementById(valueElId);
    const pctEl = document.getElementById(pctElId);

    if (valueEl) {
      const val = info.value;
      valueEl.textContent = val !== null && val !== undefined ? val.toFixed(3) : '–';
      // Color based on value (0 = red, 0.5 = grey, 1 = green)
      if (val !== null && val !== undefined) {
        const hue = val * 120; // 0 = red (0°), 0.5 = yellow-ish (60°), 1 = green (120°)
        const saturation = Math.abs(val - 0.5) * 100 + 50; // More saturated at extremes
        valueEl.style.color = `hsl(${hue}, ${saturation}%, 35%)`;
      }
    }
    if (pctEl) {
      const pct = info.percentile;
      pctEl.textContent = pct !== null && pct !== undefined ? `${pct.toFixed(0)}th percentile` : '–';
    }
  };

  formatGO('BPO', 'goBPO', 'goBPOPct');
  formatGO('CCO', 'goCCO', 'goCCOPct');
  formatGO('MFO', 'goMFO', 'goMFOPct');
}

function setText(id, text) {
  const el = document.getElementById(id);
  if (el) el.textContent = text;
}

function renderPpiSection(pair, gene1, gene2) {
  const info = pair?.ppi_network || {};
  const shared = Array.isArray(info.shared) ? info.shared : [];
  const uniqueA = Array.isArray(info.unique_gene1) ? info.unique_gene1 : [];
  const uniqueB = Array.isArray(info.unique_gene2) ? info.unique_gene2 : [];
  const geneAName = gene1.symbol || gene1.uniprot || 'Protein A';
  const geneBName = gene2.symbol || gene2.uniprot || 'Protein B';
  const toggle = document.getElementById('toggleNonShared');
  const nonSharedWrap = document.getElementById('nonSharedLists');
  const labelA = document.getElementById('ppiLabelA');
  const labelB = document.getElementById('ppiLabelB');
  if (labelA) labelA.textContent = geneAName;
  if (labelB) labelB.textContent = geneBName;

  populatePpiList('sharedPpiList', shared, MAX_SHARED_LIST, 'No shared interactors found', 'sharedPpiNote');
  populatePpiList('uniquePpiA', uniqueA, MAX_UNIQUE_LIST, `No unique partners for ${geneAName}`, 'uniquePpiANote');
  populatePpiList('uniquePpiB', uniqueB, MAX_UNIQUE_LIST, `No unique partners for ${geneBName}`, 'uniquePpiBNote');

  const hasUnique = uniqueA.length || uniqueB.length;
  if (toggle) {
    toggle.disabled = !hasUnique;
    if (!toggle.dataset.bound) {
      toggle.addEventListener('change', () => {
        showUniquePpis = !!toggle.checked;
        if (nonSharedWrap) nonSharedWrap.style.display = (showUniquePpis && hasUnique) ? '' : 'none';
        drawPpiGraph(PPI_GRAPH_DATA, showUniquePpis);
      }, { passive: true });
      toggle.dataset.bound = '1';
    }
    showUniquePpis = hasUnique ? true : false;
    toggle.checked = showUniquePpis;
  }
  if (nonSharedWrap) {
    nonSharedWrap.style.display = (hasUnique && showUniquePpis) ? '' : 'none';
  }

  if (!shared.length && !hasUnique) {
    const svg = document.getElementById('ppiNetwork');
    if (svg) svg.innerHTML = '';
    const vennSvg = document.getElementById('ppiVenn');
    if (vennSvg) vennSvg.innerHTML = '';
    const note = document.getElementById('ppiGraphNote');
    if (note) note.textContent = 'PPI data unavailable for this pair.';
    const vennStats = document.getElementById('ppiVennStats');
    if (vennStats) vennStats.style.display = 'none';
    PPI_GRAPH_DATA = null;
    setupPpiViewModeSwitch();
    return;
  }

  PPI_GRAPH_DATA = {
    gene1: geneAName,
    gene2: geneBName,
    shared,
    unique1: uniqueA,
    unique2: uniqueB
  };
  drawPpiGraph(PPI_GRAPH_DATA, showUniquePpis);
  drawPpiVenn(PPI_GRAPH_DATA);
  setupPpiViewModeSwitch();
  updatePpiView();
}

function populatePpiList(elementId, partners, limit, emptyMsg, noteId) {
  const el = document.getElementById(elementId);
  if (!el) return;
  const note = noteId ? document.getElementById(noteId) : null;
  if (!partners || !partners.length) {
    el.innerHTML = `<span class="ppi-chip empty">${emptyMsg}</span>`;
    if (note) { note.textContent = ''; note.className = 'ppi-note'; }
    return;
  }

  const renderChips = (items) => items.map((p) => {
    const label = escapeHtml(getPartnerLabel(p));
    const tooltip = escapeHtml(getPartnerTooltip(p));
    const attr = tooltip ? ` title="${tooltip}"` : '';
    return `<span class="ppi-chip"${attr}>${label}</span>`;
  }).join('');

  const hidden = partners.length - limit;
  el.innerHTML = renderChips(partners.slice(0, limit));
  el.classList.remove('expanded');

  if (note) {
    if (hidden > 0) {
      note.textContent = `+${hidden} more`;
      note.className = 'ppi-note clickable';
      note.onclick = () => {
        if (el.classList.contains('expanded')) {
          el.innerHTML = renderChips(partners.slice(0, limit));
          el.classList.remove('expanded');
          note.textContent = `+${hidden} more`;
          // Scroll back to PPI section when collapsing
          const ppiSection = document.getElementById('ppiSection');
          if (ppiSection) {
            ppiSection.scrollIntoView({ behavior: 'smooth', block: 'start' });
          }
        } else {
          el.innerHTML = renderChips(partners);
          el.classList.add('expanded');
          note.textContent = 'Show less';
        }
      };
    } else {
      note.textContent = '';
      note.className = 'ppi-note';
      note.onclick = null;
    }
  }
}

function getPartnerLabel(entry) {
  if (!entry) return 'NA';
  if (typeof entry === 'string') return entry;
  // Show only gene symbol if available, otherwise Entrez ID
  return entry.symbol || entry.id || 'NA';
}

function getPartnerShortLabel(entry) {
  if (!entry) return 'NA';
  if (typeof entry === 'string') return entry;
  // Show only gene symbol if available, otherwise Entrez ID
  return entry.symbol || entry.id || 'NA';
}

function getPartnerTooltip(entry) {
  if (!entry) return '';
  if (typeof entry === 'string') return entry;
  // Only show Entrez ID in tooltip if we have a gene symbol
  if (entry.symbol && entry.id) {
    return `${entry.symbol} (Entrez: ${entry.id})`;
  }
  return entry.id ? `Entrez: ${entry.id}` : '';
}

function escapeHtml(value) {
  if (value == null) return '';
  return String(value).replace(/[&<>"']/g, (ch) => ({
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;',
    '"': '&quot;',
    "'": '&#39;',
  }[ch]));
}

function drawPpiGraph(data, showUnique = true) {
  const svg = document.getElementById('ppiNetwork');
  const note = document.getElementById('ppiGraphNote');
  if (!svg) return;
  svg.innerHTML = '';
  if (!data) {
    if (note) note.textContent = 'PPI network unavailable for this pair.';
    return;
  }
  const width = 420;
  const height = 240;
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  const centerY = height / 2;
  const svgNS = 'http://www.w3.org/2000/svg';

  const sharedAll = Array.isArray(data.shared) ? data.shared : [];
  const unique1All = showUnique && Array.isArray(data.unique1) ? data.unique1 : [];
  const unique2All = showUnique && Array.isArray(data.unique2) ? data.unique2 : [];
  const totalPartners = sharedAll.length + unique1All.length + unique2All.length;

  // Scale partner node size based on total count
  const baseRadius = totalPartners > 100 ? 2 : totalPartners > 50 ? 3 : totalPartners > 20 ? 5 : totalPartners > 10 ? 8 : 12;
  const showLabelsAlways = totalPartners <= 15;

  const nodes = [];
  const lines = [];

  // Main gene nodes - always large with labels
  const gene1X = 100;
  const gene2X = width - 100;
  const gene1 = { key: 'gene1', label: data.gene1, x: gene1X, y: centerY, r: 20, className: 'gene-node', showLabel: true };
  const gene2 = { key: 'gene2', label: data.gene2, x: gene2X, y: centerY, r: 20, className: 'gene-node', showLabel: true };
  nodes.push(gene1, gene2);

  // Shared partners - arranged in arc between genes
  const sharedCenterX = width / 2;
  sharedAll.forEach((partner, idx) => {
    const count = sharedAll.length;
    let x, y;
    if (count <= 10) {
      // Vertical stack in center
      const frac = count === 1 ? 0.5 : idx / (count - 1);
      x = sharedCenterX;
      y = 25 + frac * (height - 50);
    } else {
      // Multiple columns for many shared partners
      const cols = Math.ceil(count / Math.ceil(height / (baseRadius * 3)));
      const col = idx % cols;
      const row = Math.floor(idx / cols);
      const rowCount = Math.ceil(count / cols);
      const colWidth = 60 / Math.max(cols - 1, 1);
      x = sharedCenterX - 30 + col * colWidth;
      y = 20 + (row / Math.max(rowCount - 1, 1)) * (height - 40);
    }
    const node = { key: `shared-${idx}`, label: getPartnerShortLabel(partner), x, y, r: baseRadius, className: 'shared-node', showLabel: showLabelsAlways };
    nodes.push(node);
    lines.push({ from: gene1, to: node, className: 'shared-link' });
    lines.push({ from: gene2, to: node, className: 'shared-link' });
  });

  // Unique partners for gene1 - fill left side with grid layout
  if (unique1All.length) {
    const count = unique1All.length;
    const availableWidth = gene1X - 25;  // Space to the left of gene1
    const availableHeight = height - 30;  // Vertical space with padding
    const nodeSize = baseRadius * 0.8;
    const spacing = Math.max(nodeSize * 2.5, 8);

    // Calculate grid dimensions based on available space and count
    const rowsInHeight = Math.floor(availableHeight / spacing);
    const neededCols = Math.ceil(count / rowsInHeight);
    const cols = Math.max(1, neededCols);
    const rowsPerCol = Math.ceil(count / cols);
    const colSpacing = availableWidth / (cols + 1);

    unique1All.forEach((partner, idx) => {
      const col = Math.floor(idx / rowsPerCol);
      const row = idx % rowsPerCol;
      const rowsInThisCol = col === cols - 1 ? count - col * rowsPerCol : rowsPerCol;

      // Position from left edge, columns go right toward gene1
      const x = 10 + colSpacing * (col + 1);
      const y = 15 + (rowsInThisCol === 1 ? availableHeight / 2 : (row / Math.max(rowsInThisCol - 1, 1)) * availableHeight);

      const node = { key: `uniqA-${idx}`, label: getPartnerShortLabel(partner), x, y, r: nodeSize, className: 'unique-node', showLabel: showLabelsAlways };
      nodes.push(node);
      lines.push({ from: gene1, to: node, className: 'non-shared' });
    });
  }

  // Unique partners for gene2 - fill right side with grid layout
  if (unique2All.length) {
    const count = unique2All.length;
    const availableWidth = width - gene2X - 25;  // Space to the right of gene2
    const availableHeight = height - 30;  // Vertical space with padding
    const nodeSize = baseRadius * 0.8;
    const spacing = Math.max(nodeSize * 2.5, 8);

    // Calculate grid dimensions based on available space and count
    const rowsInHeight = Math.floor(availableHeight / spacing);
    const neededCols = Math.ceil(count / rowsInHeight);
    const cols = Math.max(1, neededCols);
    const rowsPerCol = Math.ceil(count / cols);
    const colSpacing = availableWidth / (cols + 1);

    unique2All.forEach((partner, idx) => {
      const col = Math.floor(idx / rowsPerCol);
      const row = idx % rowsPerCol;
      const rowsInThisCol = col === cols - 1 ? count - col * rowsPerCol : rowsPerCol;

      // Position from right edge, columns go left toward gene2
      const x = width - 10 - colSpacing * (col + 1);
      const y = 15 + (rowsInThisCol === 1 ? availableHeight / 2 : (row / Math.max(rowsInThisCol - 1, 1)) * availableHeight);

      const node = { key: `uniqB-${idx}`, label: getPartnerShortLabel(partner), x, y, r: nodeSize, className: 'unique-node', showLabel: showLabelsAlways };
      nodes.push(node);
      lines.push({ from: gene2, to: node, className: 'non-shared' });
    });
  }

  // Draw lines first (behind nodes)
  lines.forEach((ln) => {
    const el = document.createElementNS(svgNS, 'line');
    el.setAttribute('x1', ln.from.x);
    el.setAttribute('y1', ln.from.y);
    el.setAttribute('x2', ln.to.x);
    el.setAttribute('y2', ln.to.y);
    el.setAttribute('class', ln.className || '');
    el.style.opacity = totalPartners > 50 ? '0.3' : totalPartners > 20 ? '0.5' : '0.7';
    svg.appendChild(el);
  });

  // Draw nodes with hover labels
  nodes.forEach((node) => {
    const group = document.createElementNS(svgNS, 'g');
    group.setAttribute('class', `ppi-node ${node.className || ''}`.trim());
    group.style.cursor = 'pointer';

    const circle = document.createElementNS(svgNS, 'circle');
    circle.setAttribute('cx', node.x);
    circle.setAttribute('cy', node.y);
    circle.setAttribute('r', node.r);
    group.appendChild(circle);

    // Label - always show for gene nodes, hover-only for partners when many
    const label = document.createElementNS(svgNS, 'text');
    const isGeneNode = node.className === 'gene-node';
    if (isGeneNode) {
      label.setAttribute('x', node.x);
      label.setAttribute('y', node.y + node.r + 14);
      label.setAttribute('text-anchor', 'middle');
      label.setAttribute('font-weight', '600');
    } else {
      label.setAttribute('x', node.x);
      label.setAttribute('y', node.y - node.r - 4);
      label.setAttribute('text-anchor', 'middle');
      label.setAttribute('font-size', '9px');
      if (!node.showLabel) {
        label.style.opacity = '0';
        label.style.transition = 'opacity 0.15s';
      }
    }
    label.textContent = node.label;
    group.appendChild(label);

    // Hover effect for partner labels
    if (!isGeneNode && !node.showLabel) {
      group.addEventListener('mouseenter', () => { label.style.opacity = '1'; });
      group.addEventListener('mouseleave', () => { label.style.opacity = '0'; });
    }

    // Add title tooltip
    const title = document.createElementNS(svgNS, 'title');
    title.textContent = node.label;
    group.appendChild(title);

    svg.appendChild(group);
  });

  if (note) {
    const counts = [];
    if (sharedAll.length) counts.push(`${sharedAll.length} shared`);
    if (unique1All.length) counts.push(`${unique1All.length} ${data.gene1}-only`);
    if (unique2All.length) counts.push(`${unique2All.length} ${data.gene2}-only`);
    if (!showUnique && (data.unique1?.length || data.unique2?.length)) {
      counts.push('unique partners hidden');
    }
    note.textContent = counts.length ? counts.join(' • ') + (totalPartners > 15 ? ' (hover for labels)' : '') : 'No PPI data available.';
  }
}

// PPI View mode switching
let currentPpiViewMode = 'venn';

function setupPpiViewModeSwitch() {
  const viewModeSelect = document.getElementById('ppiViewMode');
  const vennWrapper = document.getElementById('ppiVennWrapper');
  const networkWrapper = document.getElementById('ppiNetworkWrapper');

  if (!viewModeSelect || !vennWrapper || !networkWrapper) return;

  viewModeSelect.addEventListener('change', (e) => {
    currentPpiViewMode = e.target.value;
    updatePpiView();
  }, { passive: true });
}

function updatePpiView() {
  const vennWrapper = document.getElementById('ppiVennWrapper');
  const networkWrapper = document.getElementById('ppiNetworkWrapper');

  if (!vennWrapper || !networkWrapper) return;

  if (currentPpiViewMode === 'venn') {
    vennWrapper.style.display = 'block';
    networkWrapper.style.display = 'none';
  } else {
    vennWrapper.style.display = 'none';
    networkWrapper.style.display = 'block';
  }
}

// Hypergeometric test for PPI overlap significance
function logFactorial(n) {
  // Stirling's approximation for large n, exact for small n
  if (n < 0) return 0;
  if (n <= 1) return 0;
  if (n <= 20) {
    let result = 0;
    for (let i = 2; i <= n; i++) {
      result += Math.log(i);
    }
    return result;
  }
  // Stirling's approximation
  return n * Math.log(n) - n + 0.5 * Math.log(2 * Math.PI * n) + 1 / (12 * n);
}

function logHypergeomPMF(k, N, K, n) {
  // P(X = k) = C(K,k) * C(N-K, n-k) / C(N, n)
  // In log form: log(C(K,k)) + log(C(N-K, n-k)) - log(C(N,n))
  // log(C(a,b)) = logFactorial(a) - logFactorial(b) - logFactorial(a-b)

  const logBinom = (a, b) => {
    if (b < 0 || b > a) return -Infinity;
    return logFactorial(a) - logFactorial(b) - logFactorial(a - b);
  };

  return logBinom(K, k) + logBinom(N - K, n - k) - logBinom(N, n);
}

function hypergeomPvalue(k, N, K, n) {
  // P(X >= k) - one-tailed test for enrichment
  // k = number of shared partners (successes in sample)
  // N = total population (estimate: all unique human proteins ~20000)
  // K = total partners of gene A (successes in population)
  // n = total partners of gene B (sample size)

  let pValue = 0;
  const maxK = Math.min(K, n);

  for (let i = k; i <= maxK; i++) {
    const logP = logHypergeomPMF(i, N, K, n);
    if (logP > -700) { // Avoid underflow
      pValue += Math.exp(logP);
    }
  }

  return Math.min(pValue, 1);
}

function calculateOddsRatio(a, b, c, d) {
  // 2x2 contingency table:
  //                Interacts with B | Does not interact with B
  // Interacts with A:       a       |           b
  // Does not interact:      c       |           d
  //
  // OR = (a * d) / (b * c)

  if (b === 0 || c === 0) {
    return Infinity;
  }
  if (a === 0 || d === 0) {
    return 0;
  }
  return (a * d) / (b * c);
}

function drawPpiVenn(data) {
  const svg = document.getElementById('ppiVenn');
  const statsContainer = document.getElementById('ppiVennStats');

  if (!svg || !data) {
    if (svg) svg.innerHTML = '';
    if (statsContainer) statsContainer.style.display = 'none';
    return;
  }

  // Show stats container
  if (statsContainer) statsContainer.style.display = 'flex';

  const svgNS = 'http://www.w3.org/2000/svg';
  svg.innerHTML = '';

  const width = 420;
  const height = 280;
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);

  const shared = Array.isArray(data.shared) ? data.shared : [];
  const unique1 = Array.isArray(data.unique1) ? data.unique1 : [];
  const unique2 = Array.isArray(data.unique2) ? data.unique2 : [];

  const countA = unique1.length + shared.length;
  const countB = unique2.length + shared.length;
  const countShared = shared.length;
  const countUniqueA = unique1.length;
  const countUniqueB = unique2.length;

  // Draw Venn diagram circles with proportional sizing - bigger diagram
  const centerY = 160;

  // Scale circle radii based on partner counts (using sqrt for area proportionality)
  const maxCount = Math.max(countA, countB, 1);
  const minRadius = 45;
  const maxRadius = 100;

  // Calculate radii proportional to sqrt of count (so area is proportional to count)
  const radiusA = minRadius + (maxRadius - minRadius) * Math.sqrt(countA / maxCount);
  const radiusB = minRadius + (maxRadius - minRadius) * Math.sqrt(countB / maxCount);

  // Calculate circle positions
  let circle1X, circle2X;

  if (countShared === 0) {
    // No overlap - completely separate circles with gap
    const gap = 40;
    circle1X = width / 2 - radiusA - gap / 2;
    circle2X = width / 2 + radiusB + gap / 2;
  } else {
    // Calculate overlap distance based on shared proportion
    const minSharedCount = Math.min(countA, countB);
    const shareRatio = minSharedCount > 0 ? countShared / minSharedCount : 0;
    // overlap goes from touching (0) to significantly overlapping based on share ratio
    const overlapAmount = shareRatio * Math.min(radiusA, radiusB) * 1.2;
    const centerDistance = radiusA + radiusB - overlapAmount;
    circle1X = width / 2 - centerDistance / 2;
    circle2X = width / 2 + centerDistance / 2;
  }

  // Circle for gene A (left)
  const circle1 = document.createElementNS(svgNS, 'circle');
  circle1.setAttribute('cx', circle1X);
  circle1.setAttribute('cy', centerY);
  circle1.setAttribute('r', radiusA);
  circle1.setAttribute('fill', 'rgba(67, 160, 71, 0.3)');
  circle1.setAttribute('stroke', '#43a047');
  circle1.setAttribute('stroke-width', '2');
  svg.appendChild(circle1);

  // Circle for gene B (right)
  const circle2 = document.createElementNS(svgNS, 'circle');
  circle2.setAttribute('cx', circle2X);
  circle2.setAttribute('cy', centerY);
  circle2.setAttribute('r', radiusB);
  circle2.setAttribute('fill', 'rgba(251, 140, 0, 0.3)');
  circle2.setAttribute('stroke', '#fb8c00');
  circle2.setAttribute('stroke-width', '2');
  svg.appendChild(circle2);

  // Labels
  const addText = (x, y, text, fontSize = '14px', fontWeight = 'normal', fill = '#333') => {
    const el = document.createElementNS(svgNS, 'text');
    el.setAttribute('x', x);
    el.setAttribute('y', y);
    el.setAttribute('text-anchor', 'middle');
    el.setAttribute('font-size', fontSize);
    el.setAttribute('font-weight', fontWeight);
    el.setAttribute('fill', fill);
    el.textContent = text;
    svg.appendChild(el);
  };

  // Gene names at top with total count right below
  addText(circle1X, 28, data.gene1, '15px', '600', '#2e7d32');
  addText(circle1X, 44, `(${countA} partners)`, '11px', 'normal', '#888');
  addText(circle2X, 28, data.gene2, '15px', '600', '#e65100');
  addText(circle2X, 44, `(${countB} partners)`, '11px', 'normal', '#888');

  // Counts in circles - position based on overlap (only centered labels, no duplicates)
  if (countShared > 0) {
    // With overlap: unique counts in non-overlapping parts, shared in middle
    const overlapCenter = (circle1X + circle2X) / 2;
    const uniqueAx = circle1X - radiusA * 0.4;
    const uniqueBx = circle2X + radiusB * 0.4;
    addText(uniqueAx, centerY, countUniqueA.toString(), '20px', '700', '#2e7d32');
    addText(uniqueAx, centerY + 18, 'unique', '10px', 'normal', '#666');
    addText(overlapCenter, centerY, countShared.toString(), '20px', '700', '#5d4037');
    addText(overlapCenter, centerY + 18, 'shared', '10px', 'normal', '#666');
    addText(uniqueBx, centerY, countUniqueB.toString(), '20px', '700', '#e65100');
    addText(uniqueBx, centerY + 18, 'unique', '10px', 'normal', '#666');
  } else {
    // No overlap: counts centered in each circle, "0 shared" between
    addText(circle1X, centerY, countUniqueA.toString(), '20px', '700', '#2e7d32');
    addText(circle1X, centerY + 18, 'unique', '10px', 'normal', '#666');
    addText(width / 2, centerY, '0', '18px', '700', '#999');
    addText(width / 2, centerY + 18, 'shared', '10px', 'normal', '#999');
    addText(circle2X, centerY, countUniqueB.toString(), '20px', '700', '#e65100');
    addText(circle2X, centerY + 18, 'unique', '10px', 'normal', '#666');
  }

  // Calculate hypergeometric test
  // Estimate total human interactome size
  const TOTAL_PROTEINS = 20000; // Approximate number of human proteins

  // For hypergeometric test:
  // N = total population size
  // K = number of successes in population (gene A's partners)
  // n = number of draws (gene B's partners)
  // k = number of observed successes (shared partners)

  if (countA > 0 && countB > 0) {
    const pValue = hypergeomPvalue(countShared, TOTAL_PROTEINS, countA, countB);

    // Calculate odds ratio
    // a = shared, b = unique to A, c = unique to B, d = not interacting with either
    const a = countShared;
    const b = countUniqueA;
    const c = countUniqueB;
    const d = TOTAL_PROTEINS - countA - c; // Approximation

    const oddsRatio = calculateOddsRatio(a, b, c, d);

    // Update stats display
    const orEl = document.getElementById('vennOR');
    const pEl = document.getElementById('vennPvalue');
    const noteEl = document.getElementById('vennStatNote');

    if (orEl) {
      if (oddsRatio === Infinity) {
        orEl.textContent = '∞';
      } else if (oddsRatio === 0) {
        orEl.textContent = '0';
      } else {
        orEl.textContent = oddsRatio.toFixed(2);
      }
    }

    if (pEl) {
      if (pValue < 0.001) {
        pEl.textContent = pValue.toExponential(2);
      } else {
        pEl.textContent = pValue.toFixed(4);
      }

      // Color based on significance
      if (pValue < 0.05) {
        pEl.style.color = '#2e7d32'; // Green for significant
      } else {
        pEl.style.color = '#888';
      }
    }

    if (noteEl) {
      if (pValue < 0.001) {
        noteEl.textContent = 'Highly significant overlap (p < 0.001)';
      } else if (pValue < 0.01) {
        noteEl.textContent = 'Very significant overlap (p < 0.01)';
      } else if (pValue < 0.05) {
        noteEl.textContent = 'Significant overlap (p < 0.05)';
      } else {
        noteEl.textContent = 'Overlap not statistically significant';
      }
    }
  } else {
    // No data for test
    const orEl = document.getElementById('vennOR');
    const pEl = document.getElementById('vennPvalue');
    const noteEl = document.getElementById('vennStatNote');
    if (orEl) orEl.textContent = '–';
    if (pEl) { pEl.textContent = '–'; pEl.style.color = ''; }
    if (noteEl) noteEl.textContent = 'Insufficient data for statistical test';
  }
}

function highlightMetricSelection() {
  const cards = document.querySelectorAll('#conservationList .metric');
  cards.forEach((card) => {
    const isActive = card.dataset.metric === activeMetricKey;
    card.classList.toggle('active', isActive);
    card.setAttribute('aria-pressed', isActive ? 'true' : 'false');
  });
}

function resetMetricSelection() {
  activeMetricKey = null;
  if (boxplotChart) {
    boxplotChart.destroy();
    boxplotChart = null;
  }
  const container = document.getElementById('boxplotContainer');
  if (container) {
    container.innerHTML = `<div class="boxplot-hint">${DEFAULT_BOXPLOT_HINT}</div>`;
  }
  const details = document.getElementById('metricDetails');
  if (details) {
    details.style.display = 'none';
  }
  const title = document.getElementById('boxplotTitle');
  if (title) title.textContent = 'Select a Metric';
  const valEl = document.getElementById('detail-value');
  if (valEl) valEl.textContent = '–';
  const pctEl = document.getElementById('detail-percentile');
  if (pctEl) pctEl.textContent = '–';
  const interpEl = document.getElementById('detail-interp');
  if (interpEl) interpEl.innerHTML = '–';
  const dirEl = document.getElementById('detail-direction');
  if (dirEl) dirEl.textContent = '–';
  const fill = document.getElementById('percentile-fill');
  if (fill) fill.style.width = '0%';
  highlightMetricSelection();
}

function showBoxplotForMetric(metricKey) {
  const conservation = SUMMARY.conservation || {};
  const boxplots = SUMMARY.boxplots || {};
  
  const metricInfo = conservation[metricKey];
  const boxplotData = boxplots[metricKey];
  
  if (!metricInfo || !boxplotData) {
    document.getElementById('boxplotContainer').innerHTML = '<div class="boxplot-hint">Data not available for this metric</div>';
    return;
  }
  
  activeMetricKey = metricKey;
  highlightMetricSelection();
  
  document.getElementById('boxplotTitle').textContent = metricInfo.label || metricKey;
  document.getElementById('metricDetails').style.display = 'block';
  document.getElementById('detail-value').textContent = typeof metricInfo.value === 'number' ? metricInfo.value.toFixed(4) : '–';
  const pctVal = typeof metricInfo.percentile === 'number' ? metricInfo.percentile : null;
  document.getElementById('detail-percentile').textContent = pctVal != null ? `${pctVal.toFixed(1)}%` : '–';
  document.getElementById('percentile-fill').style.width = pctVal != null ? `${pctVal}%` : '0%';
  const dirText = metricInfo.direction_hint || (metricInfo.higher_is_more_conserved ? 'Higher values = more conserved' : 'Lower values = more conserved');
  const directionEl = document.getElementById('detail-direction');
  if (directionEl) directionEl.textContent = dirText;
  
  let interp = '';
  if (metricInfo.percentile >= 75) interp = '<span class="cons-high">Highly conserved</span> (top 25%)';
  else if (metricInfo.percentile >= 50) interp = '<span class="cons-medium">Moderately conserved</span>';
  else if (metricInfo.percentile >= 25) interp = '<span class="cons-medium">Less conserved than average</span>';
  else interp = '<span class="cons-low">Highly divergent</span> (bottom 25%)';
  document.getElementById('detail-interp').innerHTML = interp;

  drawBoxplot(metricInfo, boxplotData, metricKey);
}

function drawBoxplot(metricInfo, boxplotData, metricKey) {
  const container = document.getElementById('boxplotContainer');
  container.innerHTML = '<canvas id="boxplotCanvas" style="width:100%;height:180px;"></canvas>';

  const canvas = document.getElementById('boxplotCanvas');
  if (!canvas || typeof Chart === 'undefined') {
    container.innerHTML = '<div class="boxplot-hint">Interactive chart requires Chart.js. Metric percentiles are still listed below.</div>';
    boxplotChart = null;
    return;
  }
  const ctx = canvas.getContext('2d');
  const {q1, median, q3, whisker_low, whisker_high, pair_value} = boxplotData;

  if (boxplotChart) boxplotChart.destroy();

  // Determine if this metric should have fixed 0-1 range
  // ESM2 and ProtT5 are distance metrics that can exceed 1
  const key = metricKey || metricInfo.key || metricInfo.label || '';
  const isDistanceMetric = key.includes('esm2') || key.includes('ProtT5') ||
                           key.toLowerCase().includes('cosine');
  const useFixedRange = !isDistanceMetric;
  const xMin = useFixedRange ? 0 : whisker_low - (whisker_high - whisker_low) * 0.1;
  const xMax = useFixedRange ? 1 : whisker_high + (whisker_high - whisker_low) * 0.1;

  boxplotChart = new Chart(ctx, {
    type: 'bar',
    data: {
      labels: ['Distribution'],
      datasets: [
        { label: 'Lower', data: [q1 - whisker_low], backgroundColor: 'rgba(200,200,200,0.3)', barPercentage: 0.5 },
        { label: 'Q1-Med', data: [median - q1], backgroundColor: 'rgba(102, 126, 234, 0.4)', barPercentage: 0.5 },
        { label: 'Med-Q3', data: [q3 - median], backgroundColor: 'rgba(118, 75, 162, 0.4)', barPercentage: 0.5 },
        { label: 'Upper', data: [whisker_high - q3], backgroundColor: 'rgba(200,200,200,0.3)', barPercentage: 0.5 },
      ]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      indexAxis: 'y',
      scales: {
        x: { stacked: true, min: xMin, max: xMax, title: { display: true, text: metricInfo.label } },
        y: { stacked: true, display: false }
      },
      plugins: { legend: { display: false }, tooltip: { enabled: false } }
    },
    plugins: [{
      id: 'pairMarker',
      afterDraw: (chart) => {
        if (pair_value == null) return;
        const ctx = chart.ctx;
        const xAxis = chart.scales.x;
        const yAxis = chart.scales.y;
        const x = xAxis.getPixelForValue(pair_value);
        const y = yAxis.getPixelForValue(0);
        
        ctx.save();
        ctx.beginPath();
        ctx.arc(x, y, 10, 0, Math.PI * 2);
        ctx.fillStyle = '#ef5350';
        ctx.fill();
        ctx.strokeStyle = '#c62828';
        ctx.lineWidth = 2;
        ctx.stroke();
        ctx.fillStyle = '#333';
        ctx.font = 'bold 11px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('This pair', x, y - 16);
        ctx.font = '10px sans-serif';
        ctx.fillText(pair_value.toFixed(3), x, y + 22);
        ctx.restore();
        
        const medX = xAxis.getPixelForValue(median);
        ctx.save();
        ctx.strokeStyle = '#333';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(medX, y - 18);
        ctx.lineTo(medX, y + 18);
        ctx.stroke();
        ctx.restore();
      }
    }]
  });
}

/* ----------------- Disable unwanted external calls ----------------- *//* ----------------- Disable unwanted external calls ----------------- */
(function(){
  const blockHosts = ['molstarvolseg.ncbr.muni.cz'];
  const _fetch = window.fetch ? window.fetch.bind(window) : null;
  if (_fetch) {
    window.fetch = (input, init) => {
      try {
        const url = (typeof input === 'string') ? input : (input && input.url);
        if (url && blockHosts.some(h => url.includes(h))) {
          return Promise.resolve(new Response(JSON.stringify({ items: [] }), {
            status: 200, headers: { 'Content-Type': 'application/json' }
          }));
        }
      } catch(e) {}
      return _fetch(input, init);
    };
  }
})();

/* ----------------- Mol* wiring (main viewer) ----------------- */
let viewer = null, plugin = null, structureReady = false;
const chainIdA = 'A', chainIdB = 'B';

async function initMolstar(){
  if (viewer) return;
  const opts = {
    layoutIsExpanded:false,
    layoutShowControls:false,  // Hide Structure Tools panel by default
    layoutShowSequence:false,
    layoutShowLog:false,
    layoutShowLeftPanel:false,
    viewportShowExpand:true,   // Allow user to expand/access controls if needed
    volumeStreamingServer: ''
  };
  const v = await molstar.Viewer.create('viewer', opts);
  viewer = v;
  plugin = v.plugin;
  window.viewer = viewer;
  window.molstar = molstar;
  window.plugin = plugin;
}

async function loadPDBfromBase64(b64, resetCamera = true){
  await initMolstar();
  try { await plugin.clear(); } catch(e) {}

  const bytes = Uint8Array.from(atob(b64), c => c.charCodeAt(0));
  const blob = new Blob([bytes], {type:"chemical/x-pdb"});
  const url = URL.createObjectURL(blob);

  try {
    await viewer.loadStructureFromUrl(url, 'pdb');
    structureReady = true;
  } catch(e) {
    console.error('Failed to load structure:', e);
    structureReady = false;
  }

  if (resetCamera) {
    try { await viewer.resetCamera(); } catch(e){ plugin.canvas3d?.requestCameraReset(); }
  }
  URL.revokeObjectURL(url);
}

/* Set chain visibility in Molstar - uses component filtering approach */
async function setMolstarChainVisibility(showA, showB) {
  if (!plugin || !structureReady) return;

  console.log('setMolstarChainVisibility called: showA=', showA, 'showB=', showB);

  // If both chains should be visible, nothing to do
  if (showA && showB) {
    return;
  }

  try {
    const hierarchy = plugin.managers.structure.hierarchy.current;
    if (!hierarchy?.structures?.length) {
      console.log('No structures in hierarchy');
      return;
    }

    const struct = hierarchy.structures[0];
    const components = struct.components || [];
    console.log('Found', components.length, 'components');

    // Map unit IDs to chain IDs for the main structure
    const structData = struct.cell?.obj?.data;
    if (!structData?.units) {
      console.log('No units in structure data');
      return;
    }

    const unitIdToChain = new Map();
    for (const unit of structData.units) {
      try {
        const model = unit.model;
        if (!model?.atomicHierarchy) continue;
        const elements = unit.elements;
        if (!elements?.length) continue;
        const atomIdx = elements[0];
        const chainSegs = model.atomicHierarchy.chainAtomSegments;
        if (chainSegs?.index) {
          const chainIdx = chainSegs.index[atomIdx];
          const chainId = model.atomicHierarchy.chains.auth_asym_id.value(chainIdx);
          unitIdToChain.set(unit.id, chainId);
          console.log('Unit', unit.id, '-> chain', chainId);
        }
      } catch(e) {}
    }

    // For each component, determine which chain(s) it contains
    for (const comp of components) {
      const compData = comp.cell?.obj?.data;
      if (!compData?.units) continue;

      let compChains = new Set();
      for (const unit of compData.units) {
        const chainId = unitIdToChain.get(unit.id);
        if (chainId) compChains.add(chainId);
      }

      console.log('Component chains:', Array.from(compChains));

      // Determine visibility
      let shouldShow = true;
      if (compChains.has('A') && !compChains.has('B')) {
        shouldShow = showA;
      } else if (compChains.has('B') && !compChains.has('A')) {
        shouldShow = showB;
      }

      // Update all representations in this component
      const reps = comp.representations || [];
      for (const repr of reps) {
        if (!repr.cell?.transform?.ref) continue;
        try {
          const ref = repr.cell.transform.ref;
          const currentState = plugin.state.data.cells.get(ref)?.state;
          const newHidden = !shouldShow;
          if (currentState?.isHidden !== newHidden) {
            console.log('Setting ref', ref, 'hidden=', newHidden);
            await plugin.state.data.updateCellState(ref, { isHidden: newHidden });
          }
        } catch(e) {
          console.warn('Error updating cell state:', e);
        }
      }
    }

    console.log('Chain visibility applied: A=', showA, 'B=', showB);
  } catch(e) {
    console.warn('Could not set chain visibility:', e);
  }
}

/* ----------------- Selection + tracks + tables ----------------- */

const selection = new Map();
const trackRefs = {};
const domByUidA = {};
const domByUidB = {};

function selectionKey(chain, uid) {
  return `${chain}:${uid}`;
}

function getAllSelections() {
  return Array.from(selection.values());
}

function sanitizeRects(arr, alnLen){
  if (!Array.isArray(arr)) return [];
  const out=[];
  for (const r of arr){
    let s = Number(r.start ?? r.x ?? r.begin ?? r.from ?? 1);
    let e = Number(r.end   ?? r.to ?? r.stop  ?? r.finish ?? s);
    if (!Number.isFinite(s) || !Number.isFinite(e)) continue;
    s = Math.max(1, Math.min(alnLen, Math.floor(s)));
    e = Math.max(1, Math.min(alnLen, Math.ceil(e)));
    if (e < s) { const t = s; s = e; e = t; }
    if (e < s || e - s < 0) continue;
    const base = { x:s, start:s, begin:s };
    const color = r.color || '#999999';
    const opacity = ('opacity' in r) ? r.opacity : 1.0;
    const id = r.id;
    const label = r.label || r.name || r.type;
    out.push({ ...base, end:e, to:e, color, opacity, id, label });
  }
  return out;
}

function renderTrackSelections() {
  const sel = getAllSelections();
  const hasSel = sel.length > 0;

  const selIdsA = new Set(sel.filter(s => s.chain === chainIdA).map(s => s.id));
  const selIdsB = new Set(sel.filter(s => s.chain === chainIdB).map(s => s.id));
  const colorById = {};
  sel.forEach(s => { colorById[s.id] = s.color; });

  Object.entries(trackRefs).forEach(([name, track]) => {
    if (!track || !track._originalData) return;
    const isA = name.endsWith('A');
    const localSel = isA ? selIdsA : selIdsB;

    const newData = track._originalData.map(item => {
      const id = item.id;
      const isSelected = id && localSel.has(id);
      const baseColor = item.color || '#999999';
      return {
        ...item,
        color: isSelected ? (colorById[id] || baseColor) : baseColor,
        opacity: hasSel ? (isSelected ? 1.0 : 0.25) : 1.0
      };
    });
    track.data = newData;
  });
}

function renderTableSelections() {
  const sel = getAllSelections();
  const selKeys = new Set(sel.map(s => selectionKey(s.chain, s.id)));

  document.querySelectorAll('table#domA tbody tr, table#domB tbody tr').forEach(tr => {
    const uid = tr.getAttribute('data-uid');
    const chain = tr.getAttribute('data-chain');
    const key = selectionKey(chain, uid);
    const checked = selKeys.has(key);
    const cb = tr.querySelector('input[type="checkbox"]');
    if (cb) cb.checked = checked;
    tr.classList.toggle('selected', checked);
  });
}

let viewerLocked = true; // Default to locked
let pendingHighlightLoci = null;

function setViewerLocked(locked) {
  viewerLocked = locked;
  const btn = document.getElementById('lockViewer');
  if (btn) {
    btn.textContent = locked ? 'Unlock Hover' : 'Lock Hover';
    btn.style.background = locked ? '#ffeb3b' : '#fff';
  }
  if (!locked) pendingHighlightLoci = null;
}

function toggleViewerLock() {
  setViewerLocked(!viewerLocked);
}

// Chain visibility state
const chainVisible = { A: true, B: true };
let druggabilityFilter = 'medium+';

function getVisibleChains() {
  const chains = [];
  if (chainVisible.A) chains.push('A');
  if (chainVisible.B) chains.push('B');
  return chains;
}

async function applyChainVisibility() {
  if (!plugin) return;

  // For color modes that use special PDB variants, we need to use the appropriate PDB
  // For uniform mode, we can use chain-specific PDBs for better performance

  let pdb64;

  // If we're in a special color mode, try to use that PDB variant (contains both chains with color data)
  if (currentColorMode === 'plddt' && window.PDB64_PLDDT) {
    pdb64 = window.PDB64_PLDDT;
  } else if (currentColorMode === 'am' && window.PDB64_AM_BY_MODE && window.PDB64_AM_BY_MODE[amMode]) {
    pdb64 = window.PDB64_AM_BY_MODE[amMode];
  } else if (currentColorMode === 'aligned' && window.PDB64_ALIGNED) {
    pdb64 = window.PDB64_ALIGNED;
  } else if (currentColorMode === 'domains' && window.PDB64_DOMAINS) {
    pdb64 = window.PDB64_DOMAINS;
  }

  // Fallback: if no special variant available (or uniform mode), use chain-specific PDBs
  if (!pdb64) {
    if (chainVisible.A && chainVisible.B) {
      pdb64 = PDB64_FULL;
    } else if (chainVisible.A) {
      pdb64 = window.PDB64_A || PDB64_FULL;
    } else if (chainVisible.B) {
      pdb64 = window.PDB64_B || PDB64_FULL;
    } else {
      // Both hidden - show full anyway (will be empty view)
      pdb64 = PDB64_FULL;
    }
  }

  console.log('Applying chain visibility:', chainVisible, 'colorMode:', currentColorMode);

  // Reload with camera preservation
  await reloadViewerWith(pdb64, true);

  // Apply chain visibility via Molstar for color-mode PDBs (which contain both chains)
  if (currentColorMode && currentColorMode !== 'uniform') {
    await setMolstarChainVisibility(chainVisible.A, chainVisible.B);
    const theme = themeForColorMode(currentColorMode);
    await applyColorTheme(theme);
  }

  // Restore selections
  await renderSelections();
}

/* Apply Molstar highlighting for MAIN VIEWER - DUAL COLOR: green for A (select), pink for B (highlight) */
async function applyMolstarSelection() {
  try {
    if (!plugin || !structureReady) return;
    
    // Clear previous selections and highlights
    try {
      plugin.managers.interactivity.lociSelects.deselectAll();
    } catch(e) {}
    try {
      plugin.managers.interactivity.lociHighlights.clearHighlights();
    } catch(e) {}

    const selections = getAllSelections();
    if (selections.length === 0) {
      pendingHighlightLoci = null;
      return;
    }

    const hierarchy = plugin.managers.structure.hierarchy.current;
    if (!hierarchy?.structures?.length) return;

    const structure = hierarchy.structures[0];
    const structureData = structure.cell?.obj?.data;
    if (!structureData) return;

    const units = structureData.units || [];
    if (!units.length) return;

    console.log(`Main viewer: Structure has ${units.length} units (polymers)`);

    // Each unit represents a separate polymer
    // Unit 0 = Polymer 1 = Chain A
    // Unit 1 = Polymer 2 = Chain B
    // CRITICAL: Each unit has its own elements array, and indices in the loci must be 
    // indices INTO unit.elements, not global atom indices!
    
    const chainInfo = {}; // chainId -> { unit, unitIndex, ... }
    
    for (let unitIdx = 0; unitIdx < units.length; unitIdx++) {
      const unit = units[unitIdx];
      try {
        const model = unit.model;
        if (!model) continue;
        
        const chains = model.atomicHierarchy?.chains;
        const chainAtomSegments = model.atomicHierarchy?.chainAtomSegments;
        const residueAtomSegments = model.atomicHierarchy?.residueAtomSegments;
        
        if (!chains || !chainAtomSegments || !residueAtomSegments) continue;
        
        // Get the elements (atom indices) that this unit contains
        const unitElements = unit.elements;
        if (!unitElements || unitElements.length === 0) continue;
        
        // Find which chain the first atom of this unit belongs to
        const firstAtomIdx = unitElements[0];
        const chainOffsets = chainAtomSegments.offsets;
        
        let unitChainId = null;
        for (let ci = 0; ci < chainOffsets.length - 1; ci++) {
          if (firstAtomIdx >= chainOffsets[ci] && firstAtomIdx < chainOffsets[ci + 1]) {
            unitChainId = chains.label_asym_id.value(ci);
            break;
          }
        }
        
        if (!unitChainId) {
          console.warn(`Unit ${unitIdx}: could not determine chain ID`);
          continue;
        }
        
        // Build a mapping from protein residue number (1-based) to local atom indices
        const resOffsets = residueAtomSegments.offsets;
        const unitResidueToLocalAtoms = new Map();
        
        // Determine the first residue in the chain (to compute 1-based protein residue)
        let minGlobalResIdx = Infinity;
        
        // First pass: find all global residue indices in this unit
        for (let localIdx = 0; localIdx < unitElements.length; localIdx++) {
          const globalAtomIdx = unitElements[localIdx];
          for (let ri = 0; ri < resOffsets.length - 1; ri++) {
            if (globalAtomIdx >= resOffsets[ri] && globalAtomIdx < resOffsets[ri + 1]) {
              if (ri < minGlobalResIdx) minGlobalResIdx = ri;
              break;
            }
          }
        }
        
        // Second pass: for each atom, map to protein residue and collect local indices
        for (let localIdx = 0; localIdx < unitElements.length; localIdx++) {
          const globalAtomIdx = unitElements[localIdx];
          for (let ri = 0; ri < resOffsets.length - 1; ri++) {
            if (globalAtomIdx >= resOffsets[ri] && globalAtomIdx < resOffsets[ri + 1]) {
              const proteinRes = ri - minGlobalResIdx + 1;
              if (!unitResidueToLocalAtoms.has(proteinRes)) {
                unitResidueToLocalAtoms.set(proteinRes, []);
              }
              unitResidueToLocalAtoms.get(proteinRes).push(localIdx);
              break;
            }
          }
        }
        
        const residueCount = unitResidueToLocalAtoms.size;
        console.log(`Unit ${unitIdx}: chain ${unitChainId}, ${unitElements.length} atoms, ${residueCount} residues`);
        
        chainInfo[unitChainId] = {
          unit: unit,
          unitIndex: unitIdx,
          residueCount: residueCount,
          unitResidueToLocalAtoms: unitResidueToLocalAtoms
        };
        
      } catch(e) {
        console.warn(`Error analyzing unit ${unitIdx}:`, e);
      }
    }

    // Fallback if detection failed
    if (Object.keys(chainInfo).length === 0) {
      console.warn('Chain detection failed, using simple unit index fallback');
      for (let unitIdx = 0; unitIdx < Math.min(2, units.length); unitIdx++) {
        const chainId = unitIdx === 0 ? 'A' : 'B';
        const unit = units[unitIdx];
        const unitElements = unit.elements;
        if (!unitElements) continue;
        
        const unitResidueToLocalAtoms = new Map();
        const atomsPerResidue = 8;
        let proteinRes = 1;
        for (let i = 0; i < unitElements.length; i += atomsPerResidue) {
          const localAtoms = [];
          for (let j = i; j < Math.min(i + atomsPerResidue, unitElements.length); j++) {
            localAtoms.push(j);
          }
          unitResidueToLocalAtoms.set(proteinRes, localAtoms);
          proteinRes++;
        }
        
        chainInfo[chainId] = {
          unit: unit,
          unitIndex: unitIdx,
          residueCount: unitResidueToLocalAtoms.size,
          unitResidueToLocalAtoms: unitResidueToLocalAtoms
        };
      }
    }

    // Separate selections by chain for dual-color highlighting
    const elementsChainA = [];
    const elementsChainB = [];
    
    for (const sel of selections) {
      try {
        const targetChain = sel.chain; // 'A' or 'B'
        const info = chainInfo[targetChain];
        
        if (!info) {
          console.warn(`No chain info found for chain ${targetChain}`);
          continue;
        }

        const { unit, unitIndex, unitResidueToLocalAtoms } = info;
        
        const startProteinRes = sel.start; // 1-based
        const endProteinRes = sel.end;     // 1-based
        
        console.log(`Selection "${sel.name}" chain ${targetChain}: protein res ${startProteinRes}-${endProteinRes}`);
        
        const localAtomIndices = [];
        
        for (let proteinRes = startProteinRes; proteinRes <= endProteinRes; proteinRes++) {
          const localAtoms = unitResidueToLocalAtoms.get(proteinRes);
          if (localAtoms) {
            localAtomIndices.push(...localAtoms);
          }
        }

        if (localAtomIndices.length > 0) {
          console.log(`  -> ${localAtomIndices.length} LOCAL atoms selected for unit ${unitIndex} (chain ${targetChain})`);
          
          localAtomIndices.sort((a, b) => a - b);
          
          if (targetChain === 'A') {
            elementsChainA.push({ unit, indices: localAtomIndices });
          } else {
            elementsChainB.push({ unit, indices: localAtomIndices });
          }
        } else {
          console.warn(`  -> No atoms found for selection`);
        }
      } catch (e) {
        console.warn(`Failed to process selection ${sel.name}:`, e);
      }
    }
    
    // Apply chain A selections using lociSelects (GREEN color)
    if (elementsChainA.length > 0) {
      const lociA = {
        kind: 'element-loci',
        structure: structureData,
        elements: elementsChainA
      };
      try {
        plugin.managers.interactivity.lociSelects.select({ loci: lociA });
        console.log(`Applied SELECT (green) for chain A with ${elementsChainA.length} element groups`);
      } catch (e) {
        console.warn('Failed to apply chain A selection:', e);
      }
    }
    
    // Apply chain B selections using lociHighlights (PINK color)
    if (elementsChainB.length > 0) {
      const lociB = {
        kind: 'element-loci',
        structure: structureData,
        elements: elementsChainB
      };
      
      pendingHighlightLoci = lociB;
      
      try {
        plugin.managers.interactivity.lociHighlights.highlight({ loci: lociB });
        console.log(`Applied HIGHLIGHT (pink) for chain B with ${elementsChainB.length} element groups`);
      } catch (e) {
        console.warn('Failed to apply chain B highlight:', e);
      }
    } else {
      pendingHighlightLoci = null;
    }
  } catch (e) {
    console.error('applyMolstarSelection failed:', e);
  }
}

/* Re-apply the chain B highlight (called when hover would normally clear it) */
function reapplyChainBHighlight() {
  if (pendingHighlightLoci && viewerLocked) {
    try {
      plugin.managers.interactivity.lociHighlights.highlight({ loci: pendingHighlightLoci });
    } catch(e) {}
  }
}

/* Setup hover interception to maintain highlights when locked */
let hoverInterceptionSetup = false;
let hoverReapplyInterval = null;

function setupHoverInterception() {
  if (hoverInterceptionSetup) {
    console.log('Hover interception already setup, skipping');
    return;
  }

  if (!plugin?.canvas3d?.interaction?.hover) {
    console.warn('Cannot setup hover interception - interaction.hover not available');
    return;
  }

  try {
    // Subscribe to hover events
    plugin.canvas3d.interaction.hover.subscribe((e) => {
      if (viewerLocked && pendingHighlightLoci) {
        // Immediately reapply highlight
        reapplyChainBHighlight();
      }
    });

    // Also set up continuous reapplication when locked
    // This ensures highlights persist even during continuous mouse movement
    if (!hoverReapplyInterval) {
      hoverReapplyInterval = setInterval(() => {
        if (viewerLocked && pendingHighlightLoci) {
          reapplyChainBHighlight();
        }
      }, 50); // Reapply every 50ms when locked
    }

    hoverInterceptionSetup = true;
    console.log('Hover interception setup complete with continuous reapplication');
  } catch(e) {
    console.warn('Failed to setup hover interception:', e);
  }
}

async function initializeHighlightColors() {
  try {
    if (plugin?.canvas3d?.setProps) {
      await plugin.canvas3d.setProps({
        marking: {
          selectColor: { r: 0.26, g: 0.63, b: 0.28 },
          highlightColor: { r: 0.91, g: 0.12, b: 0.39 }
        }
      });
      return true;
    }
  } catch(e) {}
  return false;
}

window.toggleViewerLock = toggleViewerLock;
window.setViewerLocked = setViewerLocked;

/* =============================================================================
   PDBe VIEWER - v6 with simplified chain selection
   ============================================================================= */

let pdbeViewer = null;
let pdbePlugin = null;
let pdbeStructureReady = false;
let currentPdbeEntry = null;
let currentPdbeIndex = -1;

async function initPdbeMolstar() {
  if (pdbeViewer) return;
  const container = document.getElementById('pdbeViewer');
  if (!container) {
    console.error('PDBe viewer container not found');
    return;
  }
  
  try {
    const v = await molstar.Viewer.create('pdbeViewer', {
      layoutIsExpanded: false,
      layoutShowControls: false,  // Hide Structure Tools panel by default
      layoutShowSequence: false,
      layoutShowLog: false,
      layoutShowLeftPanel: false,
      viewportShowExpand: true,   // Allow user to expand if needed
      volumeStreamingServer: ''
    });
    pdbeViewer = v;
    pdbePlugin = v.plugin;
    
    try {
      await pdbePlugin.canvas3d.setProps({
        marking: {
          selectColor: { r: 0.26, g: 0.63, b: 0.28 },
          highlightColor: { r: 0.91, g: 0.12, b: 0.39 }
        }
      });
    } catch(e) {
      console.warn('Could not set PDBe highlight colors:', e);
    }
    
    console.log('PDBe Molstar viewer initialized');
  } catch(e) {
    console.error('Failed to initialize PDBe viewer:', e);
  }
}

async function applyGreyColoring() {
  if (!pdbePlugin || !pdbeStructureReady) return;
  
  try {
    const structures = pdbePlugin.managers.structure.hierarchy.current.structures;
    if (!structures || structures.length === 0) return;
    
    for (const struct of structures) {
      const components = struct.components || [];
      for (const comp of components) {
        if (comp.representations) {
          for (const repr of comp.representations) {
            try {
              const update = pdbePlugin.state.data.build().to(repr.cell)
                .update(old => {
                  return {
                    ...old,
                    colorTheme: {
                      name: 'uniform',
                      params: { value: 0xaaaaaa }
                    }
                  };
                });
              await update.commit();
            } catch(e) {}
          }
        }
      }
    }
    console.log('Applied grey coloring to structure');
  } catch(e) {
    console.warn('Could not apply grey coloring:', e);
  }
}

async function loadPdbeStructureFromBase64(b64, format = 'pdb') {
  await initPdbeMolstar();
  if (!pdbeViewer || !pdbePlugin) {
    console.error('PDBe viewer not available');
    return false;
  }

  try { await pdbePlugin.clear(); } catch (e) {}

  if (!b64 || b64.length < 10) {
    console.warn('No valid base64 data provided');
    pdbeStructureReady = false;
    return false;
  }

  try {
    const bytes = Uint8Array.from(atob(b64), c => c.charCodeAt(0));
    const mimeType = format === 'mmcif' || format === 'cif' ? 'chemical/x-mmcif' : 'chemical/x-pdb';
    const blob = new Blob([bytes], { type: mimeType });
    const url = URL.createObjectURL(blob);

    const molstarFormat = (format === 'mmcif' || format === 'cif') ? 'mmcif' : 'pdb';
    
    await pdbeViewer.loadStructureFromUrl(url, molstarFormat);
    pdbeStructureReady = true;
    
    URL.revokeObjectURL(url);
    
    try { 
      await pdbeViewer.resetCamera(); 
    } catch(e) { 
      pdbePlugin.canvas3d?.requestCameraReset(); 
    }
    
    setTimeout(() => applyGreyColoring(), 500);
    
    console.log('PDBe structure loaded successfully');
    return true;
  } catch (e) {
    console.error('Failed to load PDBe structure:', e);
    pdbeStructureReady = false;
    return false;
  }
}

async function fetchAndLoadPdbeStructure(pdbId) {
  await initPdbeMolstar();
  if (!pdbeViewer || !pdbePlugin || !pdbId) return false;

  try { await pdbePlugin.clear(); } catch(e) {}

  const pid = pdbId.toLowerCase();
  
  const sources = [
    { url: `https://files.rcsb.org/download/${pid}.pdb`, format: 'pdb' },
    { url: `https://www.ebi.ac.uk/pdbe/entry-files/download/pdb${pid}.ent`, format: 'pdb' },
    { url: `https://files.rcsb.org/download/${pid}.cif`, format: 'mmcif' },
  ];

  for (const src of sources) {
    try {
      console.log(`Trying to fetch from: ${src.url}`);
      await pdbeViewer.loadStructureFromUrl(src.url, src.format);
      pdbeStructureReady = true;
      console.log(`Successfully loaded ${pdbId} from ${src.url}`);
      
      try { await pdbeViewer.resetCamera(); } catch(e) {}
      setTimeout(() => applyGreyColoring(), 500);
      
      return true;
    } catch(e) {
      console.warn(`Failed to load from ${src.url}:`, e.message);
    }
  }
  
  console.error(`Could not load structure ${pdbId} from any source`);
  pdbeStructureReady = false;
  return false;
}

function getTargetChainIds(entry) {
  if (!entry) return new Set();
  
  const chains = entry.chains || entry.chain_id || entry.chainId || '';
  const targetChains = new Set();
  
  if (Array.isArray(chains)) {
    chains.forEach(c => {
      if (c && typeof c === 'string') {
        targetChains.add(c.trim());
      }
    });
  } else if (typeof chains === 'string') {
    chains.split(/[,\s]+/).forEach(c => {
      const trimmed = c.trim();
      if (trimmed) targetChains.add(trimmed);
    });
  }
  
  console.log('Target chain IDs from PDBe entry:', Array.from(targetChains));
  return targetChains;
}

/**
 * Highlight chains by auth_asym_id - simplified version that trusts PDBe data
 * Iterates all units and checks each one's chain assignment
 */
async function highlightProteinChains(entry) {
  if (!pdbePlugin || !pdbeStructureReady || !entry) {
    console.warn('PDBe viewer not ready for highlighting');
    return;
  }

  const targetChains = getTargetChainIds(entry);
  if (targetChains.size === 0) {
    console.warn('No target chains specified');
    updatePdbeLegend('No chains specified for highlighting.');
    return;
  }

  try {
    try {
      pdbePlugin.managers.interactivity.lociSelects.deselectAll();
      pdbePlugin.managers.interactivity.lociHighlights.clearHighlights();
    } catch(e) {}

    const structures = pdbePlugin.managers.structure.hierarchy.current.structures;
    if (!structures || structures.length === 0) {
      console.warn('No structures loaded');
      return;
    }

    const struct = structures[0];
    const structureData = struct.cell?.obj?.data;
    if (!structureData) {
      console.warn('No structure data');
      return;
    }

    // Log structure info for debugging
    console.log('Structure units:', structureData.units?.length);
    
    const units = structureData.units || [];
    const matchingElements = [];
    const matchedChains = new Set();
    const allChainsFound = new Set();

    for (let unitIdx = 0; unitIdx < units.length; unitIdx++) {
      const unit = units[unitIdx];
      try {
        const model = unit.model;
        if (!model) continue;

        const chains = model.atomicHierarchy?.chains;
        const chainAtomSegments = model.atomicHierarchy?.chainAtomSegments;
        
        if (!chains || !chainAtomSegments) continue;

        const unitElements = unit.elements;
        if (!unitElements || !unitElements.length) continue;

        // Get auth_asym_id for this unit by checking its first atom
        const firstAtomIdx = unitElements[0];
        const chainOffsets = chainAtomSegments.offsets;
        
        let unitAuthChainId = null;
        
        for (let ci = 0; ci < chainOffsets.length - 1; ci++) {
          if (firstAtomIdx >= chainOffsets[ci] && firstAtomIdx < chainOffsets[ci + 1]) {
            unitAuthChainId = chains.auth_asym_id?.value(ci);
            allChainsFound.add(unitAuthChainId);
            break;
          }
        }
        
        if (!unitAuthChainId) continue;
        
        // Check if this unit's chain matches our targets
        if (targetChains.has(unitAuthChainId)) {
          console.log(`Unit ${unitIdx}: chain ${unitAuthChainId} matches target`);
          matchedChains.add(unitAuthChainId);
          
          // Add all atoms in this unit
          const indices = [];
          for (let i = 0; i < unitElements.length; i++) {
            indices.push(i);
          }
          matchingElements.push({ unit, indices });
        }
      } catch(e) {
        console.warn(`Error analyzing unit ${unitIdx}:`, e);
      }
    }

    console.log('All chains found in structure:', Array.from(allChainsFound));
    console.log('Matched chains:', Array.from(matchedChains));

    if (matchingElements.length > 0) {
      const loci = {
        kind: 'element-loci',
        structure: structureData,
        elements: matchingElements
      };

      pdbePlugin.managers.interactivity.lociSelects.select({ loci });
      
      console.log(`Highlighted ${matchingElements.length} units for chains: ${Array.from(matchedChains).join(', ')}`);
      updatePdbeLegend(`<strong>Highlighted:</strong> Chain(s) <strong>${Array.from(matchedChains).join(', ')}</strong> shown in <span style="color:#43a047;font-weight:bold">green</span>.`);
    } else {
      console.warn('No matching chains found. Available chains:', Array.from(allChainsFound));
      updatePdbeLegend(`<strong>Note:</strong> Could not find chains ${Array.from(targetChains).join(', ')}. Available: ${Array.from(allChainsFound).join(', ')}`);
    }

  } catch(e) {
    console.error('Error highlighting protein chains:', e);
  }
}

/**
 * Highlight ligands - finds non-polymer residues (HETATM records)
 * Uses Mol* selection queries for more reliable detection
 */
async function highlightLigands() {
  if (!pdbePlugin || !pdbeStructureReady) {
    console.warn('PDBe viewer not ready for highlighting');
    return;
  }

  try {
    // Clear existing selections/highlights
    try {
      pdbePlugin.managers.interactivity.lociSelects.deselectAll();
      pdbePlugin.managers.interactivity.lociHighlights.clearHighlights();
    } catch(e) {}

    const structures = pdbePlugin.managers.structure.hierarchy.current.structures;
    if (!structures || structures.length === 0) {
      console.warn('No structures loaded');
      return;
    }

    const struct = structures[0];
    const structureData = struct.cell?.obj?.data;
    if (!structureData) {
      console.warn('No structure data');
      return;
    }

    // Standard amino acids and common non-ligands to exclude
    const excludeResidues = new Set([
      'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
      'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL',
      'MSE','SEC','PYL',  // Modified amino acids
      'A','C','G','U','T','DA','DC','DG','DT','DU',  // Nucleotides
      'HOH','WAT','H2O','DOD','SOL'  // Water
    ]);

    const units = structureData.units || [];
    const matchingElements = [];
    const ligandNames = new Set();

    for (let unitIdx = 0; unitIdx < units.length; unitIdx++) {
      const unit = units[unitIdx];
      try {
        const model = unit.model;
        if (!model) continue;

        // Access residue data
        const residueNames = model.atomicHierarchy?.residues?.label_comp_id;
        const atomSegments = model.atomicHierarchy?.residueAtomSegments;
        
        if (!residueNames || !atomSegments) continue;

        const unitElements = unit.elements;
        if (!unitElements || unitElements.length === 0) continue;

        const ligandAtomIndices = [];
        const resOffsets = atomSegments.offsets;

        // Check each atom in the unit
        for (let i = 0; i < unitElements.length; i++) {
          const atomIdx = unitElements[i];
          
          // Find which residue this atom belongs to
          for (let ri = 0; ri < resOffsets.length - 1; ri++) {
            if (atomIdx >= resOffsets[ri] && atomIdx < resOffsets[ri + 1]) {
              const resName = residueNames.value(ri);
              const resNameUpper = resName ? resName.toUpperCase() : '';
              
              // Check if this is a ligand (not standard AA, nucleotide, or water)
              if (resNameUpper && !excludeResidues.has(resNameUpper)) {
                ligandAtomIndices.push(i);
                ligandNames.add(resNameUpper);
              }
              break;
            }
          }
        }

        if (ligandAtomIndices.length > 0) {
          matchingElements.push({ unit, indices: ligandAtomIndices });
        }

      } catch(e) {
        console.warn(`Error analyzing unit ${unitIdx} for ligands:`, e);
      }
    }

    console.log('Ligand detection result:', {
      unitsChecked: units.length,
      matchingElements: matchingElements.length,
      ligandNames: Array.from(ligandNames)
    });

    if (matchingElements.length > 0 && ligandNames.size > 0) {
      const loci = {
        kind: 'element-loci',
        structure: structureData,
        elements: matchingElements
      };

      // Use select instead of highlight for more visible coloring
      pdbePlugin.managers.interactivity.lociSelects.select({ loci });
      
      const ligandList = Array.from(ligandNames).join(', ');
      console.log(`Highlighted ligands: ${ligandList}`);
      updatePdbeLegend(`<strong>Highlighted:</strong> Ligands (<strong>${ligandList}</strong>) shown in <span style="color:#43a047;font-weight:bold">green</span>.`);
    } else {
      console.warn('No ligands found in structure');
      updatePdbeLegend(`<strong>Note:</strong> No ligands found in this structure.`);
    }

  } catch(e) {
    console.error('Error highlighting ligands:', e);
  }
}

async function clearPdbeHighlights() {
  if (!pdbePlugin) return;
  
  try {
    pdbePlugin.managers.interactivity.lociSelects.deselectAll();
    pdbePlugin.managers.interactivity.lociHighlights.clearHighlights();
    
    updatePdbeLegend(`Complex shown in <strong>grey</strong>. Click buttons to highlight protein of interest (<span style="color:#43a047">green</span>) or ligands (<span style="color:#e91e63">pink</span>).`);
  } catch(e) {
    console.warn('Error clearing highlights:', e);
  }
}

function updatePdbeLegend(html) {
  const legend = document.getElementById('pdbeLegend');
  if (legend) legend.innerHTML = html;
}

async function syncCameraOrientationFromMainViewer() {
  if (!plugin || !pdbePlugin || !structureReady || !pdbeStructureReady) {
    console.warn('Both viewers must be ready to sync camera');
    return;
  }
  
  try {
    const mainCamera = plugin.canvas3d?.camera;
    const pdbeCamera = pdbePlugin.canvas3d?.camera;
    
    if (!mainCamera || !pdbeCamera) {
      console.warn('Camera not available');
      return;
    }
    
    const mainSnapshot = mainCamera.getSnapshot();
    const pdbeSnapshot = pdbeCamera.getSnapshot();
    
    // Calculate main camera's view direction (normalized)
    const mainDir = [
      mainSnapshot.target[0] - mainSnapshot.position[0],
      mainSnapshot.target[1] - mainSnapshot.position[1],
      mainSnapshot.target[2] - mainSnapshot.position[2]
    ];
    const mainDist = Math.sqrt(mainDir[0]**2 + mainDir[1]**2 + mainDir[2]**2);
    if (mainDist > 0) {
      mainDir[0] /= mainDist;
      mainDir[1] /= mainDist;
      mainDir[2] /= mainDist;
    }
    
    // Calculate current PDBe camera distance from target
    const pdbeDir = [
      pdbeSnapshot.target[0] - pdbeSnapshot.position[0],
      pdbeSnapshot.target[1] - pdbeSnapshot.position[1],
      pdbeSnapshot.target[2] - pdbeSnapshot.position[2]
    ];
    const pdbeDist = Math.sqrt(pdbeDir[0]**2 + pdbeDir[1]**2 + pdbeDir[2]**2);
    
    // Apply main camera's direction to PDBe camera, keeping PDBe's distance and target
    const newState = {
      ...pdbeSnapshot,
      up: mainSnapshot.up,
      position: [
        pdbeSnapshot.target[0] - mainDir[0] * pdbeDist,
        pdbeSnapshot.target[1] - mainDir[1] * pdbeDist,
        pdbeSnapshot.target[2] - mainDir[2] * pdbeDist
      ]
    };
    
    await pdbeCamera.setState(newState, 300);
    
    console.log('Camera orientation synced from main viewer');
  } catch(e) {
    console.warn('Could not sync camera orientation:', e);
  }
}

function updatePdbeInfoBox(entry) {
  if (!entry) {
    document.getElementById('pdbePdbId').textContent = '-';
    document.getElementById('pdbeSourceAcc').textContent = '-';
    document.getElementById('pdbeChains').textContent = '-';
    document.getElementById('pdbeCoverage').textContent = '-';
    document.getElementById('pdbeResolution').textContent = '-';
    document.getElementById('pdbeLigands').textContent = '-';
    return;
  }

  const pdbId = (entry.pdb_id || entry.pdbId || '').toUpperCase();
  const sourceAcc = entry.source_acc || entry.sourceAcc || entry.uniprot_acc || '';
  const chains = entry.chains || entry.chain_id || entry.chainId || '';
  const coverage = entry.coverage;
  const resolution = entry.resolution;
  const ligandSummary = entry.ligandSummary || entry.ligand_summary || '';

  document.getElementById('pdbePdbId').textContent = pdbId || '-';
  document.getElementById('pdbeSourceAcc').textContent = sourceAcc || '-';
  document.getElementById('pdbeChains').textContent = 
    Array.isArray(chains) ? chains.join(', ') : (chains || '-');
  document.getElementById('pdbeCoverage').textContent = 
    (typeof coverage === 'number') ? (coverage * 100).toFixed(1) + '%' : '-';
  document.getElementById('pdbeResolution').textContent = 
    (typeof resolution === 'number') ? resolution.toFixed(2) + ' Å' : '-';
  document.getElementById('pdbeLigands').textContent = ligandSummary || 'None reported';
}

async function loadPdbeByIndex(idx) {
  console.log('loadPdbeByIndex called with:', idx, 'PDBe_COMPLEXES length:', PDBe_COMPLEXES.length);
  if (typeof idx === 'string') idx = parseInt(idx, 10);
  if (!Number.isFinite(idx) || idx < 0 || idx >= PDBe_COMPLEXES.length) {
    console.warn('Invalid PDBe index:', idx);
    return;
  }

  const entry = PDBe_COMPLEXES[idx];
  console.log('Loading PDBe entry:', entry);
  if (!entry) return;

  currentPdbeEntry = entry;
  currentPdbeIndex = idx;

  // Update button text with gene name
  const sourceAcc = entry.source_acc || entry.sourceAcc || '';
  const geneName = getGeneNameForAccession(sourceAcc);
  const proteinNameSpan = document.getElementById('pdbeProteinName');
  if (proteinNameSpan) {
    proteinNameSpan.textContent = geneName || 'Protein';
  }

  updatePdbeInfoBox(entry);
  updatePdbeLegend(`Loading structure...`);

  const b64 = entry.coord_b64 || entry.coordB64 || entry.pdb_b64 || entry.pdbB64 || '';
  const format = entry.coord_format || entry.coordFormat || 'pdb';
  const pdbId = (entry.pdb_id || entry.pdbId || '').toLowerCase();

  let loaded = false;

  if (b64 && b64.length > 100) {
    console.log(`Loading PDBe structure ${pdbId} from base64 (format: ${format})`);
    loaded = await loadPdbeStructureFromBase64(b64, format);
  }

  if (!loaded && pdbId) {
    console.log(`Fetching PDBe structure ${pdbId} from remote`);
    loaded = await fetchAndLoadPdbeStructure(pdbId);
  }

  if (loaded) {
    updatePdbeLegend(`Complex shown in <strong>grey</strong>. <strong>${geneName || 'Protein'}</strong> highlighted in <span style="color:#43a047">green</span>.`);
    // Auto-highlight the protein of interest
    setTimeout(async () => {
      await highlightProteinChains(entry);
      updateHighlightButtonState(true);
    }, 300);
  } else {
    console.error('Failed to load structure for entry:', entry);
    updatePdbeLegend(`<strong>Error:</strong> Failed to load structure.`);
  }
}

// Helper to get gene name for a UniProt accession
function getGeneNameForAccession(acc) {
  if (!SUMMARY || !acc) return null;
  const acc1 = SUMMARY.gene1?.uniprot;
  const acc2 = SUMMARY.gene2?.uniprot;
  if (acc === acc1) return SUMMARY.gene1?.symbol;
  if (acc === acc2) return SUMMARY.gene2?.symbol;
  return null;
}

function updateHighlightButtonState(highlighted) {
  isProteinHighlighted = highlighted;
  const btn = document.getElementById('pdbeHighlightProtein');
  if (btn) {
    if (highlighted) {
      btn.style.backgroundColor = '#43a047';
      btn.style.color = 'white';
      btn.style.borderColor = '#388e3c';
    } else {
      btn.style.backgroundColor = '';
      btn.style.color = '';
      btn.style.borderColor = '';
    }
  }
}

function setupPdbeControls() {
  const structSelect = document.getElementById('pdbeStructSelect');
  const highlightProteinBtn = document.getElementById('pdbeHighlightProtein');
  const clearBtn = document.getElementById('pdbeClearHighlight');
  const syncBtn = document.getElementById('pdbeSyncCamera');
  const centerBtn = document.getElementById('pdbeCenter');

  if (!structSelect) return;

  structSelect.innerHTML = '';

  if (!PDBe_COMPLEXES || PDBe_COMPLEXES.length === 0) {
    const opt = document.createElement('option');
    opt.value = '';
    opt.textContent = 'No experimental structures available';
    structSelect.appendChild(opt);
    structSelect.disabled = true;
    if (highlightProteinBtn) highlightProteinBtn.disabled = true;
    updatePdbeInfoBox(null);
    return;
  }

  structSelect.disabled = false;
  if (highlightProteinBtn) highlightProteinBtn.disabled = false;

  PDBe_COMPLEXES.forEach((entry, idx) => {
    const pdbId = (entry.pdb_id || entry.pdbId || '').toUpperCase();
    const sourceAcc = entry.source_acc || entry.sourceAcc || '';
    const chains = entry.chains || entry.chain_id || '';
    const resolution = entry.resolution;
    const ligandSummary = entry.ligandSummary || entry.ligand_summary || '';

    let label = pdbId;
    if (sourceAcc) label += ` [${sourceAcc}]`;
    if (chains) {
      const chainStr = Array.isArray(chains) ? chains.join(',') : chains;
      label += ` ch:${chainStr}`;
    }
    if (typeof resolution === 'number') {
      label += ` ${resolution.toFixed(1)}Å`;
    }
    if (ligandSummary) {
      const shortLig = ligandSummary.length > 15 ? ligandSummary.slice(0, 15) + '…' : ligandSummary;
      label += ` (${shortLig})`;
    }

    const opt = document.createElement('option');
    opt.value = String(idx);
    opt.textContent = label;
    structSelect.appendChild(opt);
  });

  structSelect.addEventListener('change', async (e) => {
    await loadPdbeByIndex(e.target.value);
  }, { passive: true });

  if (highlightProteinBtn) {
    highlightProteinBtn.addEventListener('click', async () => {
      if (currentPdbeEntry) {
        if (isProteinHighlighted) {
          await clearPdbeHighlights();
          updateHighlightButtonState(false);
        } else {
          await highlightProteinChains(currentPdbeEntry);
          updateHighlightButtonState(true);
        }
      }
    }, { passive: true });
  }

  if (clearBtn) {
    clearBtn.addEventListener('click', async () => {
      await clearPdbeHighlights();
      updateHighlightButtonState(false);
    }, { passive: true });
  }

  if (syncBtn) {
    syncBtn.addEventListener('click', async () => {
      await syncCameraOrientationFromMainViewer();
    }, { passive: true });
  }

  if (centerBtn) {
    centerBtn.addEventListener('click', () => {
      if (pdbeViewer && pdbeStructureReady) {
        try {
          pdbeViewer.resetCamera();
        } catch(e) {
          pdbePlugin?.canvas3d?.requestCameraReset();
        }
      }
    }, { passive: true });
  }

  if (structSelect.options.length > 0) {
    structSelect.selectedIndex = 0;
    // Load first structure after a small delay to ensure viewer is ready
    setTimeout(async () => {
      await loadPdbeByIndex(0);
    }, 500);
  }
}

async function renderSelections() {
  renderTrackSelections();
  renderTableSelections();
  await applyMolstarSelection();
}

function toggleFeature(dom, chain) {
  if (!dom || !dom.uid) return;
  const key = selectionKey(chain, dom.uid);
  if (selection.has(key)) {
    selection.delete(key);
  } else {
    const color = (dom.type === 'Cavity') ? '#ff7d45' : '#ffdb13';
    selection.set(key, {
      id: dom.uid,
      chain,
      start: parseInt(dom.start, 10),
      end: parseInt(dom.end, 10),
      color,
      name: dom.label || dom.name || dom.type || 'region'
    });
  }
  renderSelections();
}

function applyAmMode(mode){
  const modes = AM_MODES;
  if (!modes.includes(mode)) mode = modes[0] || 'raw';
  amMode = mode;

  const alnLen = Math.max(1, (DATA.qaln || '').length);

  if (amTrackA) {
    const rectsA = (DATA.amAlignedARectsByMode && DATA.amAlignedARectsByMode[mode]) || [];
    amTrackA.data = sanitizeRects(rectsA, alnLen);
  }
  if (amTrackB) {
    const rectsB = (DATA.amAlignedBRectsByMode && DATA.amAlignedBRectsByMode[mode]) || [];
    amTrackB.data = sanitizeRects(rectsB, alnLen);
  }

  if (damTrack) {
    const dr = (DATA.damAlignedRectsByMode && DATA.damAlignedRectsByMode[mode]) || [];
    damTrack.data = sanitizeRects(dr, alnLen);
  }

  const matAByMode = DATA.amMatrixA_rectsByMode || {};
  const matBByMode = DATA.amMatrixB_rectsByMode || {};
  const matAForMode = matAByMode[mode] || {};
  const matBForMode = matBByMode[mode] || {};

  amMatrixTracksA.forEach(({track, aa}) => {
    const rects = matAForMode[aa] || [];
    track.data = sanitizeRects(rects, alnLen);
  });
  amMatrixTracksB.forEach(({track, aa}) => {
    const rects = matBForMode[aa] || [];
    track.data = sanitizeRects(rects, alnLen);
  });
}

function addRow(tbl, label, el, h){
  const tr=document.createElement('tr');
  const td1=document.createElement('td'); td1.className='rowlbl'; td1.textContent=label;
  const td2=document.createElement('td'); td2.append(el);
  if (h) el.setAttribute('height', String(h));
  tr.append(td1,td2); tbl.append(tr);
  return { row: tr, labelCell: td1, trackCell: td2 };
}

function buildSeq(){
  const mgr = document.getElementById('mgr'); mgr.innerHTML='';
  const tbl = document.createElement('table'); tbl.className='rtab'; mgr.append(tbl);

  const alnLen = Math.max(1, (DATA.qaln||'').length);

  amMatrixTracksA = [];
  amMatrixTracksB = [];

  const nav = document.createElement('nightingale-navigation');
  addRow(tbl, '', nav, 40);

  const seqA = document.createElement('nightingale-sequence');
  addRow(tbl, DATA.g1+' ('+DATA.a1+') aligned', seqA, 28);

  const domA = document.createElement('nightingale-track');
  addRow(tbl, 'Domains '+DATA.g1, domA, 16); trackRefs['domA'] = domA;

  const disorderA = document.createElement('nightingale-track');
  addRow(tbl, 'Disordered '+DATA.g1, disorderA, 16); trackRefs['disorderA'] = disorderA;

  const tedA = document.createElement('nightingale-track');
  addRow(tbl, 'TED '+DATA.g1, tedA, 16); trackRefs['tedA'] = tedA;

  const cavA = document.createElement('nightingale-track');
  addRow(tbl, 'Cavities '+DATA.g1, cavA, 16); trackRefs['cavA'] = cavA;

  const amA = document.createElement('nightingale-track');
  const amARow = addRow(tbl, 'AlphaMissense '+DATA.g1, amA, 18);
  amARow.row.classList.add('am-main-row');
  amARow.labelCell.classList.add('am-main-label');
  amTrackA = amA;

  const amMatrixRowsA = [];
  const toggleA = document.createElement('button');
  toggleA.textContent = '▼';
  toggleA.className = 'am-matrix-toggle';
  amARow.labelCell.appendChild(toggleA);

  AA_ORDER.forEach(aa => {
    const trk = document.createElement('nightingale-track');
    const rowObj = addRow(tbl, aa, trk, 8);
    rowObj.row.style.display = 'none';
    amMatrixRowsA.push(rowObj.row);
    amMatrixTracksA.push({ track: trk, aa });
  });

  const parentA = amARow.row.parentElement;
  amMatrixRowsA.forEach(r => {
    parentA.insertBefore(r, amARow.row);
  });

  toggleA.addEventListener('click', () => {
    const hidden = amMatrixRowsA.length && amMatrixRowsA[0].style.display === 'none';
    amMatrixRowsA.forEach(r => { r.style.display = hidden ? '' : 'none'; });
    toggleA.textContent = hidden ? '▲' : '▼';
  }, { passive:true });

  const dam = document.createElement('nightingale-track');
  const damRow = addRow(tbl, 'ΔAM', dam, 18);
  damRow.row.classList.add('am-main-row');
  damRow.labelCell.classList.add('am-main-label');
  damTrack = dam;

  const amB = document.createElement('nightingale-track');
  const amBRow = addRow(tbl, 'AlphaMissense '+DATA.g2, amB, 18);
  amBRow.row.classList.add('am-main-row');
  amBRow.labelCell.classList.add('am-main-label');
  amTrackB = amB;

  const amMatrixRowsB = [];
  const toggleB = document.createElement('button');
  toggleB.textContent = '▼';
  toggleB.className = 'am-matrix-toggle';
  amBRow.labelCell.appendChild(toggleB);

  AA_ORDER.forEach(aa => {
    const trk = document.createElement('nightingale-track');
    const rowObj = addRow(tbl, aa, trk, 8);
    rowObj.row.style.display = 'none';
    amMatrixRowsB.push(rowObj.row);
    amMatrixTracksB.push({ track: trk, aa });
  });

  toggleB.addEventListener('click', () => {
    const hidden = amMatrixRowsB.length && amMatrixRowsB[0].style.display === 'none';
    amMatrixRowsB.forEach(r => { r.style.display = hidden ? '' : 'none'; });
    toggleB.textContent = hidden ? '▲' : '▼';
  }, { passive:true });

  const seqB = document.createElement('nightingale-sequence');
  addRow(tbl, DATA.g2+' ('+DATA.a2+') aligned', seqB, 28);

  const domB = document.createElement('nightingale-track');
  addRow(tbl, 'Domains '+DATA.g2, domB, 16); trackRefs['domB'] = domB;

  const disorderB = document.createElement('nightingale-track');
  addRow(tbl, 'Disordered '+DATA.g2, disorderB, 16); trackRefs['disorderB'] = disorderB;

  const tedB = document.createElement('nightingale-track');
  addRow(tbl, 'TED '+DATA.g2, tedB, 16); trackRefs['tedB'] = tedB;

  const cavB = document.createElement('nightingale-track');
  addRow(tbl, 'Cavities '+DATA.g2, cavB, 16); trackRefs['cavB'] = cavB;

  requestAnimationFrame(()=>{
    const allTracks = [
      nav, seqA,
      domA, disorderA, tedA, cavA,
      amA, dam, amB,
      seqB,
      domB, disorderB, tedB, cavB
    ];
    amMatrixTracksA.forEach(obj => allTracks.push(obj.track));
    amMatrixTracksB.forEach(obj => allTracks.push(obj.track));

    allTracks.forEach(el => {
      if (!el) return;
      el.length = alnLen;
      el.displaystart = 1;
      el.displayend = alnLen;
    });

    seqA.data = DATA.qaln || '';
    seqB.data = DATA.taln || '';

    amA.setAttribute('shape','rectangle');
    amB.setAttribute('shape','rectangle');
    dam.setAttribute('shape','rectangle');

    [domA, disorderA, tedA, cavA, domB, disorderB, tedB, cavB].forEach(track => {
      track.setAttribute('shape','roundRectangle');
      track.setAttribute('show-label','');
    });

    domA.data       = sanitizeRects(DATA.domA_alnRects||[], alnLen);          domA._originalData = [...domA.data];
    disorderA.data  = sanitizeRects(DATA.disorderA_alnRects||[], alnLen);     disorderA._originalData = [...disorderA.data];
    tedA.data       = sanitizeRects(DATA.tedA_alnRects||[], alnLen);          tedA._originalData = [...tedA.data];
    cavA.data       = sanitizeRects(DATA.cavA_alnRects||[], alnLen);          cavA._originalData = [...cavA.data];

    domB.data       = sanitizeRects(DATA.domB_alnRects||[], alnLen);          domB._originalData = [...domB.data];
    disorderB.data  = sanitizeRects(DATA.disorderB_alnRects||[], alnLen);     disorderB._originalData = [...disorderB.data];
    tedB.data       = sanitizeRects(DATA.tedB_alnRects||[], alnLen);          tedB._originalData = [...tedB.data];
    cavB.data       = sanitizeRects(DATA.cavB_alnRects||[], alnLen);          cavB._originalData = [...cavB.data];

    applyAmMode(amMode);
    applyCavityFilter(); // Apply default druggability filter (medium+)
  });

  function attachDomainClick(track, chain) {
    track.addEventListener('click', (ev)=>{
      const feat = ev?.detail?.feature;
      if (!feat || !feat.id) return;
      const domMap = (chain === chainIdA) ? domByUidA : domByUidB;
      const dom = domMap[feat.id];
      if (!dom) return;
      toggleFeature(dom, chain);
    }, {passive:true});
  }

  attachDomainClick(domA, chainIdA);
  attachDomainClick(disorderA, chainIdA);
  attachDomainClick(tedA, chainIdA);
  attachDomainClick(cavA, chainIdA);
  attachDomainClick(domB, chainIdB);
  attachDomainClick(disorderB, chainIdB);
  attachDomainClick(tedB, chainIdB);
  attachDomainClick(cavB, chainIdB);
}

function setupPdbeCollapse(){
  const btn = document.getElementById('pdbeCollapseBtn');
  const body = document.getElementById('pdbeCardBody');
  const card = document.getElementById('pdbeCard');
  if (!btn || !body) return;
  btn.addEventListener('click', () => {
    const collapsed = body.classList.toggle('collapsed');
    btn.setAttribute('aria-expanded', (!collapsed).toString());
    btn.textContent = collapsed ? 'Expand' : 'Collapse';
    if (card) card.classList.toggle('is-collapsed', collapsed);
  }, {passive:true});
}

// Generic collapsible section setup
function setupCollapsibleSection(btnId, bodyId, sectionId) {
  const btn = document.getElementById(btnId);
  const body = document.getElementById(bodyId);
  const section = sectionId ? document.getElementById(sectionId) : null;
  if (!btn || !body) return;
  btn.addEventListener('click', () => {
    const collapsed = body.classList.toggle('collapsed');
    btn.setAttribute('aria-expanded', (!collapsed).toString());
    btn.textContent = collapsed ? 'Expand' : 'Collapse';
    if (section) section.classList.toggle('is-collapsed', collapsed);
  }, {passive:true});
}

// Hide section if it has no data
function hideSection(sectionId) {
  const section = document.getElementById(sectionId);
  if (section) {
    section.classList.add('section-hidden');
  }
}

// Show section
function showSection(sectionId) {
  const section = document.getElementById(sectionId);
  if (section) {
    section.classList.remove('section-hidden');
  }
}

// Setup all collapsible sections
function setupAllCollapsibleSections() {
  // Family nav
  setupCollapsibleSection('familyNavCollapseBtn', 'familyNavBody', 'familyNav');
  // Family features
  setupCollapsibleSection('familyFeaturesCollapseBtn', 'familyFeaturesBody', 'familyFeaturesSection');
  // Summary
  setupCollapsibleSection('summaryCollapseBtn', 'summaryBody', 'summarySection');
  // Similarity search
  setupCollapsibleSection('simSearchCollapseBtn', 'simSearchBody', 'similaritySearchSection');
  // Structure
  setupCollapsibleSection('structureCollapseBtn', 'structureBody', 'structureSection');
  // Alignment
  setupCollapsibleSection('alignmentCollapseBtn', 'alignmentBody', 'alignmentSection');
  // Domains
  setupCollapsibleSection('domainsCollapseBtn', 'domainsBody', 'domainsSection');
  // Domain pairs
  setupCollapsibleSection('domainPairsCollapseBtn', 'domainPairsBody', 'domainPairsSection');
}

// Check data availability and hide empty sections
function updateSectionVisibility() {
  // PDBe section - hide if no experimental structures
  if (!PDBe_COMPLEXES || PDBe_COMPLEXES.length === 0) {
    hideSection('pdbeCard');
  } else {
    showSection('pdbeCard');
  }

  // Structure section - hide if no PDB data
  if (!PDB64_FULL) {
    hideSection('structureSection');
    hideSection('alignmentSection');
  } else {
    showSection('structureSection');
    showSection('alignmentSection');
  }

  // Domains section - hide if no domains
  const hasDomainsA = DATA && DATA.domainsA && DATA.domainsA.length > 0;
  const hasDomainsB = DATA && DATA.domainsB && DATA.domainsB.length > 0;
  if (!hasDomainsA && !hasDomainsB) {
    hideSection('domainsSection');
  } else {
    showSection('domainsSection');
  }

  // Domain pairs section - hide if no domain pairs
  const hasDomainPairs = DATA && DATA.domPairs && DATA.domPairs.length > 0;
  if (!hasDomainPairs) {
    hideSection('domainPairsSection');
  } else {
    showSection('domainPairsSection');
  }

  // Family features section - hide if no family features data
  const hasFamilyFeatures = SUMMARY && SUMMARY.family_features && Object.keys(SUMMARY.family_features).length > 0;
  if (!hasFamilyFeatures) {
    hideSection('familyFeaturesSection');
  } else {
    showSection('familyFeaturesSection');
  }

  // Similarity search section - hide if no similarity search data
  const hasSimilaritySearch = SUMMARY && SUMMARY.similarity_search && Object.keys(SUMMARY.similarity_search).length > 0;
  if (!hasSimilaritySearch) {
    hideSection('similaritySearchSection');
  } else {
    showSection('similaritySearchSection');
  }

  // Update sidebar navigation to match
  updateSidebarVisibility();
}

// Update sidebar links visibility based on section visibility
function updateSidebarVisibility() {
  const sidebarLinks = document.querySelectorAll('.sidebar-nav a[data-section]');
  sidebarLinks.forEach(link => {
    const sectionId = link.dataset.section;
    const section = document.getElementById(sectionId);
    const listItem = link.closest('li');
    if (section && listItem) {
      if (section.classList.contains('section-hidden') ||
          (section.style.display === 'none' && !section.classList.contains('section-hidden'))) {
        listItem.style.display = 'none';
      } else {
        listItem.style.display = '';
      }
    }
  });
}

function shouldShowDomain(d) {
  // Apply druggability filter only to cavities
  if (d.type !== 'Cavity' && d.raw_type !== 'Cavity') return true;

  const dg = (d.druggability || '').toLowerCase();

  if (druggabilityFilter === 'all') return true;
  if (druggabilityFilter === 'high') return dg === 'high';
  if (druggabilityFilter === 'medium+') return dg === 'high' || dg === 'medium';

  return true;
}

function shouldShowCavityRect(rect) {
  // Filter cavity rectangles based on druggability
  if (!rect.druggability) return true; // No druggability info, show it

  const dg = (rect.druggability || '').toLowerCase();

  if (druggabilityFilter === 'all') return true;
  if (druggabilityFilter === 'high') return dg === 'high';
  if (druggabilityFilter === 'medium+') return dg === 'high' || dg === 'medium';

  return true;
}

function applyCavityFilter() {
  // Apply druggability filter to Nightingale cavity tracks
  const cavA = trackRefs['cavA'];
  const cavB = trackRefs['cavB'];

  if (cavA && cavA._originalData) {
    cavA.data = cavA._originalData.filter(shouldShowCavityRect);
  }

  if (cavB && cavB._originalData) {
    cavB.data = cavB._originalData.filter(shouldShowCavityRect);
  }

  console.log('Applied cavity filter:', druggabilityFilter);
}

/**
 * Compute simple global sequence alignment (Needleman-Wunsch)
 * Used as fallback when pre-computed alignment is not available
 */
function computeSimpleAlignment(seq1, seq2) {
  const GAP = -2;
  const MATCH = 2;
  const MISMATCH = -1;

  const m = seq1.length;
  const n = seq2.length;

  // Initialize score matrix
  const score = Array(m + 1).fill(null).map(() => Array(n + 1).fill(0));
  for (let i = 0; i <= m; i++) score[i][0] = i * GAP;
  for (let j = 0; j <= n; j++) score[0][j] = j * GAP;

  // Fill score matrix
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const match = score[i-1][j-1] + (seq1[i-1] === seq2[j-1] ? MATCH : MISMATCH);
      const del = score[i-1][j] + GAP;
      const ins = score[i][j-1] + GAP;
      score[i][j] = Math.max(match, del, ins);
    }
  }

  // Traceback
  let qaln = '', taln = '';
  let i = m, j = n;
  const aligned_cols = [];
  let col = 0;

  while (i > 0 || j > 0) {
    if (i > 0 && j > 0 && score[i][j] === score[i-1][j-1] + (seq1[i-1] === seq2[j-1] ? MATCH : MISMATCH)) {
      qaln = seq1[i-1] + qaln;
      taln = seq2[j-1] + taln;
      aligned_cols.unshift([col, i-1, j-1]);
      i--; j--;
    } else if (i > 0 && score[i][j] === score[i-1][j] + GAP) {
      qaln = seq1[i-1] + qaln;
      taln = '-' + taln;
      aligned_cols.unshift([col, i-1, null]);
      i--;
    } else {
      qaln = '-' + qaln;
      taln = seq2[j-1] + taln;
      aligned_cols.unshift([col, null, j-1]);
      j--;
    }
    col++;
  }

  // Reindex columns
  aligned_cols.forEach((c, idx) => c[0] = idx);

  // Calculate identity
  let matches = 0;
  for (let k = 0; k < qaln.length; k++) {
    if (qaln[k] !== '-' && taln[k] !== '-' && qaln[k] === taln[k]) matches++;
  }
  const identity = matches / Math.min(m, n);

  return { qaln, taln, aligned_cols, identity };
}

async function switchAlignmentMethod(method) {
  if (!DATA || !DATA.PAIR) return;

  const btn = document.getElementById('alignmentMethod');
  if (btn) btn.disabled = true;

  try {
    if (method === 'sequence') {
      console.log('Switching to sequence alignment...');

      // Use pre-computed sequence alignment from report data, or compute client-side
      let seqData = DATA.seqAlignment;
      if (!seqData && DATA.seq1 && DATA.seq2) {
        // Compute simple sequence alignment client-side
        seqData = computeSimpleAlignment(DATA.seq1, DATA.seq2);
      }
      if (!seqData) {
        throw new Error('Sequence alignment not available');
      }
      console.log(`Sequence alignment: ${seqData.aligned_cols.length} columns, identity=${(seqData.identity * 100).toFixed(1)}%`);

      // Store original structural alignment
      if (!window.STRUCTURAL_ALIGNMENT) {
        window.STRUCTURAL_ALIGNMENT = {
          qaln: DATA.qaln,
          taln: DATA.taln,
          qposByCol: DATA.qposByCol,
          tposByCol: DATA.tposByCol,
          aligned_cols: [], // would need to reconstruct if needed
        };
      }

      // Update DATA with sequence alignment
      DATA.qaln = seqData.qaln;
      DATA.taln = seqData.taln;

      // Build new position maps
      const qposByCol = {};
      const tposByCol = {};
      for (const [col, qpos, tpos] of seqData.aligned_cols) {
        if (qpos !== null) qposByCol[col] = qpos;
        if (tpos !== null) tposByCol[col] = tpos;
      }
      DATA.qposByCol = qposByCol;
      DATA.tposByCol = tposByCol;

      // Reload all tracks with new alignment
      await reloadTracksWithAlignment();

      console.log('Switched to sequence alignment');

    } else if (method === 'structural') {
      console.log('Switching to structural alignment...');

      if (window.STRUCTURAL_ALIGNMENT) {
        // Restore structural alignment
        DATA.qaln = window.STRUCTURAL_ALIGNMENT.qaln;
        DATA.taln = window.STRUCTURAL_ALIGNMENT.taln;
        DATA.qposByCol = window.STRUCTURAL_ALIGNMENT.qposByCol;
        DATA.tposByCol = window.STRUCTURAL_ALIGNMENT.tposByCol;

        // Reload all tracks
        await reloadTracksWithAlignment();

        console.log('Switched back to structural alignment');
      } else {
        // Already on structural, just reload
        window.location.reload();
      }
    }

  } catch (e) {
    console.error('Failed to switch alignment:', e);
    alert(`Failed to switch alignment: ${e.message}`);
    if (btn) btn.value = method === 'sequence' ? 'structural' : 'sequence';
  } finally {
    if (btn) btn.disabled = false;
  }
}

async function reloadTracksWithAlignment() {
  // Reload just the Nightingale section with new alignment
  const seqwrap = document.getElementById('seqwrap');
  if (!seqwrap) return;

  // Show loading overlay on the section
  const loadingDiv = document.createElement('div');
  loadingDiv.style.cssText = 'position:fixed;inset:0;display:flex;align-items:center;justify-content:center;background:rgba(255,255,255,0.9);z-index:10000;font-size:18px;color:#333;font-weight:600;';
  loadingDiv.innerHTML = '<div>Recalculating tracks with new alignment...<div style="font-size:14px;margin-top:8px;font-weight:normal;color:#666;">This may take a moment</div></div>';
  document.body.appendChild(loadingDiv);

  try {
    // Recalculate all alignment-based rectangles
    await recalculateAlignmentTracks();

    // Clear and rebuild Nightingale manager
    seqwrap.innerHTML = '<nightingale-manager id="mgr" class="mgr"></nightingale-manager>';

    // Small delay to let DOM update
    await new Promise(resolve => setTimeout(resolve, 100));

    // Rebuild the Nightingale tracks
    buildSeq();

    console.log('Tracks reloaded with new alignment');
  } catch (e) {
    console.error('Failed to reload tracks:', e);
    throw e;
  } finally {
    if (loadingDiv.parentNode) {
      loadingDiv.parentNode.removeChild(loadingDiv);
    }
  }
}

async function recalculateAlignmentTracks() {
  // Recalculate all alignment-based tracks with new position mapping
  const alnLen = DATA.qaln ? DATA.qaln.length : 0;

  console.log(`Recalculating tracks for alignment length: ${alnLen}`);

  // Build position-to-column maps from qaln/taln
  let qpos = 0;
  let tpos = 0;
  const qposToCol = {};
  const tposToCol = {};

  for (let col = 0; col < alnLen; col++) {
    const qaa = DATA.qaln[col];
    const taa = DATA.taln[col];

    if (qaa !== '-') {
      qpos++;
      qposToCol[qpos] = col + 1; // 1-indexed columns
    }

    if (taa !== '-') {
      tpos++;
      tposToCol[tpos] = col + 1;
    }
  }

  console.log(`Position maps built: qpos=${qpos}, tpos=${tpos}`);

  // Helper to remap domains/features to alignment columns
  function remapFeaturesToAlignment(features, posToCol, defaultColor) {
    const rects = [];

    for (const feat of features) {
      const start = feat.start || 0;
      const end = feat.end || 0;

      // Find all columns that overlap this feature
      const cols = [];
      for (let pos = start; pos <= end; pos++) {
        if (posToCol[pos]) {
          cols.push(posToCol[pos]);
        }
      }

      if (cols.length > 0) {
        // Create contiguous segments
        cols.sort((a, b) => a - b);
        let segStart = cols[0];
        let prev = cols[0];

        for (let i = 1; i <= cols.length; i++) {
          const curr = i < cols.length ? cols[i] : null;

          if (curr === null || curr !== prev + 1) {
            // End of segment
            rects.push({
              start: segStart,
              end: prev,
              color: feat.color || defaultColor,
              label: feat.label || '',
              id: feat.uid || `${start}-${end}`,
              druggability: feat.druggability
            });

            if (curr !== null) {
              segStart = curr;
            }
          }

          if (curr !== null) {
            prev = curr;
          }
        }
      }
    }

    return rects;
  }

  // Recalculate domain/disorder/cavity rectangles
  if (DATA.domainsA) {
    const doms = DATA.domainsA.filter(d => d.type !== 'CAV' && d.type !== 'Cavity');
    const cavs = DATA.domainsA.filter(d => d.type === 'CAV' || d.type === 'Cavity');

    DATA.domA_alnRects = remapFeaturesToAlignment(doms, qposToCol, '#2ca02c');
    DATA.cavA_alnRects = remapFeaturesToAlignment(cavs, qposToCol, '#ff7d45');
  }

  if (DATA.domainsB) {
    const doms = DATA.domainsB.filter(d => d.type !== 'CAV' && d.type !== 'Cavity');
    const cavs = DATA.domainsB.filter(d => d.type === 'CAV' || d.type === 'Cavity');

    DATA.domB_alnRects = remapFeaturesToAlignment(doms, tposToCol, '#2ca02c');
    DATA.cavB_alnRects = remapFeaturesToAlignment(cavs, tposToCol, '#ff7d45');
  }

  // Recalculate AM tracks if we have bfactors data
  if (DATA.bfactorsA && DATA.bfactorsB && DATA.amModes) {
    console.log('Recalculating AM tracks...');

    // Helper to remap per-residue scores to alignment axis
    function remapScoresToAlignment(scores, posToCol, alnLen) {
      const alignedScores = new Array(alnLen).fill(null);
      for (let pos = 1; pos <= scores.length; pos++) {
        const col = posToCol[pos];
        if (col !== undefined && col >= 1 && col <= alnLen) {
          alignedScores[col - 1] = scores[pos - 1];
        }
      }
      return alignedScores;
    }

    // Remap bfactors to alignment axis
    const amA_aligned = remapScoresToAlignment(DATA.bfactorsA, qposToCol, alnLen);
    const amB_aligned = remapScoresToAlignment(DATA.bfactorsB, tposToCol, alnLen);

    // Helper to create rect segments from score array
    function scoreArrayToRects(scoreArr, colorFn) {
      const rects = [];
      let curColor = null;
      let segStart = null;

      for (let i = 0; i < scoreArr.length; i++) {
        const score = scoreArr[i];
        const col = i + 1;

        let color = null;
        if (score !== null && score !== undefined) {
          color = colorFn(score);
        }

        if (color === null) {
          // Close current segment if any
          if (curColor !== null && segStart !== null) {
            rects.push({start: segStart, end: col - 1, color: curColor});
          }
          curColor = null;
          segStart = null;
        } else if (curColor === null) {
          // Start new segment
          curColor = color;
          segStart = col;
        } else if (color !== curColor) {
          // Color changed - close current, start new
          rects.push({start: segStart, end: col - 1, color: curColor});
          curColor = color;
          segStart = col;
        }
      }

      // Close final segment
      if (curColor !== null && segStart !== null) {
        rects.push({start: segStart, end: alnLen, color: curColor});
      }

      return rects;
    }

    // AM color bands (same as Python)
    function bandColorAM(v) {
      if (v === null || v === undefined) return null;
      if (v < 0.2) return '#dddddd';  // Benign - grey
      if (v < 0.4) return '#bbbbbb';  // Low
      if (v < 0.7) return '#ff7d45';  // Medium - orange
      return '#d62728';                // High - red
    }

    // Delta-AM color
    function damColor(v) {
      if (v === null || v === undefined) return null;
      if (v <= 0.10) return '#dddddd';
      if (v <= 0.30) return '#bbbbbb';
      if (v <= 0.50) return '#ffdb13';  // Yellow
      if (v <= 0.70) return '#ff7d45';  // Orange
      return '#d62728';                 // Red
    }

    // Normalize scores based on mode
    function makeNormalizer(values, mode) {
      const validVals = values.filter(v => typeof v === 'number' && !isNaN(v));
      if (validVals.length === 0 || mode === 'raw') {
        return v => (typeof v === 'number' && !isNaN(v)) ? v : null;
      }

      if (mode === 'minmax') {
        const min = Math.min(...validVals);
        const max = Math.max(...validVals);
        if (max <= min) return v => 0.5;
        return v => (typeof v === 'number' && !isNaN(v)) ? (v - min) / (max - min) : null;
      }

      if (mode === 'percentile') {
        const sorted = [...validVals].sort((a, b) => a - b);
        return v => {
          if (typeof v !== 'number' || isNaN(v)) return null;
          let idx = 0;
          while (idx < sorted.length && sorted[idx] <= v) idx++;
          return sorted.length > 1 ? idx / sorted.length : 0.5;
        };
      }

      if (mode === 'zscore') {
        const mean = validVals.reduce((a, b) => a + b, 0) / validVals.length;
        const variance = validVals.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / validVals.length;
        const stddev = Math.sqrt(variance) || 1e-8;
        return v => {
          if (typeof v !== 'number' || isNaN(v)) return null;
          const z = (v - mean) / stddev;
          return 1.0 / (1.0 + Math.exp(-z));  // Sigmoid
        };
      }

      return v => (typeof v === 'number' && !isNaN(v)) ? v : null;
    }

    // Recalculate for each mode
    DATA.amAlignedARectsByMode = {};
    DATA.amAlignedBRectsByMode = {};
    DATA.damAlignedRectsByMode = {};

    for (const mode of DATA.amModes) {
      const normA = makeNormalizer(DATA.bfactorsA, mode);
      const normB = makeNormalizer(DATA.bfactorsB, mode);

      // Normalize aligned scores
      const normAmA = amA_aligned.map(v => normA(v));
      const normAmB = amB_aligned.map(v => normB(v));

      // Generate rectangles
      DATA.amAlignedARectsByMode[mode] = scoreArrayToRects(normAmA, bandColorAM);
      DATA.amAlignedBRectsByMode[mode] = scoreArrayToRects(normAmB, bandColorAM);

      // Delta-AM (only where both aligned)
      const damArr = new Array(alnLen).fill(null);
      for (let i = 0; i < alnLen; i++) {
        const vA = normAmA[i];
        const vB = normAmB[i];
        if (typeof vA === 'number' && typeof vB === 'number' && !isNaN(vA) && !isNaN(vB)) {
          damArr[i] = Math.abs(vA - vB);
        }
      }
      DATA.damAlignedRectsByMode[mode] = scoreArrayToRects(damArr, damColor);
    }

    console.log('AM tracks recalculated for new alignment');
  } else {
    console.warn('AM tracks not recalculated - bfactors data not available');
  }

  console.log('Recalculated alignment tracks');
}

function fillDomainTables(){
  const tA = document.querySelector('#domA tbody'); tA.innerHTML='';
  const tB = document.querySelector('#domB tbody'); tB.innerHTML='';

  const addRowDom = (tb, d, chain) => {
    const tr=document.createElement('tr');
    tr.className='clickable';
    tr.setAttribute('data-uid', d.uid || '');
    tr.setAttribute('data-chain', chain);

    const tdSel = document.createElement('td');
    const cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.addEventListener('change', (ev)=>{
      ev.stopPropagation();
      toggleFeature(d, chain);
    }, {passive:true});
    tdSel.append(cb);

    const tdName = document.createElement('td');
    tdName.textContent = d.label || d.name || d.type || '';

    const tdRange = document.createElement('td');
    tdRange.textContent = `${d.start}-${d.end}`;

    const tdCav = document.createElement('td');
    if (d.type === 'Cavity' || d.raw_type === 'Cavity') {
      const ds = d.drug_score || d.drugscore || '';
      const dg = d.druggability || '';
      tdCav.textContent = (ds || dg) ? `${ds || ''} ${dg || ''}`.trim() : '';
    } else {
      tdCav.textContent = '';
    }

    tr.append(tdSel, tdName, tdRange, tdCav);
    tr.addEventListener('click', ()=>{
      toggleFeature(d, chain);
    }, {passive:true});
    tb.append(tr);
  };

  (DATA.domainsA||[]).filter(shouldShowDomain).forEach(d=>addRowDom(tA,d,chainIdA));
  (DATA.domainsB||[]).filter(shouldShowDomain).forEach(d=>addRowDom(tB,d,chainIdB));

  renderTableSelections();
}

async function fillDomPairs(){
  const tb = document.querySelector('#domPairs tbody'); tb.innerHTML='';
  if (!Array.isArray(DATA.domPairs) || !DATA.domPairs.length){
    tb.innerHTML = '<tr><td colspan="5" class="small">No domain sub-alignments computed.</td></tr>';
    return;
  }
  DATA.domPairs.forEach((r)=>{
    const tr=document.createElement('tr'); tr.className='clickable';
    tr.innerHTML = `<td>${r.Aname} ${r.Arng}</td><td>${r.Bname} ${r.Brng}</td><td>${r.fident!=null?r.fident.toFixed(1)+'%':'–'}</td><td>${r.tm!=null?r.tm.toFixed(3):'–'}</td><td>${r.damPct!=null?r.damPct.toFixed(1)+'%':'–'}</td>`;
    tr.addEventListener('click', async ()=>{
      await reloadViewerWith(r.pdb64);
      document.getElementById('tmScore').textContent = (r.tm!=null ? r.tm.toFixed(3) : '–');
      const title = `Domain: ${r.Aname} ${r.Arng} × ${r.Bname} ${r.Brng}`;
      document.getElementById('contextTitle').textContent = title;

      selection.clear();
      const domA = (DATA.domainsA||[]).find(d => (d.label===r.Aname) || (`${d.start}-${d.end}`===r.Arng));
      const domB = (DATA.domainsB||[]).find(d => (d.label===r.Bname) || (`${d.start}-${d.end}`===r.Brng));
      if (domA) selection.set(selectionKey(chainIdA, domA.uid), {
        id: domA.uid, chain: chainIdA,
        start: parseInt(domA.start,10), end: parseInt(domA.end,10),
        color: '#ffdb13', name: domA.label||domA.name
      });
      if (domB) selection.set(selectionKey(chainIdB, domB.uid), {
        id: domB.uid, chain: chainIdB,
        start: parseInt(domB.start,10), end: parseInt(domB.end,10),
        color: '#ff7d45', name: domB.label||domB.name
      });
      await renderSelections();
    }, {passive:true});
    tb.append(tr);
  });
}

let currentColorMode = 'uniform';

async function applyColorTheme(theme){
  if (!plugin || !structureReady) return;
  const hierarchy = plugin.managers.structure.hierarchy.current;
  if (!hierarchy || !hierarchy.structures || !hierarchy.structures.length) return;

  const update = plugin.state.data.build();
  for (const struct of hierarchy.structures) {
    const components = struct.components || [];
    for (const comp of components) {
      const reps = comp.representations || [];
      for (const repr of reps) {
        update.to(repr.cell).update(old => ({
          ...old,
          colorTheme: theme
        }));
      }
    }
  }
  await update.commit();
}

function themeForColorMode(mode){
  const m = (mode || 'uniform').toLowerCase();
  if (m === 'plddt') {
    // pLDDT coloring using Molstar's official thresholds and colors
    // Server encodes B-factors as bin indices (0-3):
    //   0 = ≤50: Very Low confidence
    //   1 = 51-70: Low confidence
    //   2 = 71-90: Confident
    //   3 = >90: Very High confidence
    // NOTE: Molstar's uncertainty theme has `reverse: true`, so colors are applied
    // in REVERSE order: high B-factor → first color, low B-factor → last color
    // So we reverse the color array: [blue, cyan, yellow, orange]
    return {
      name: 'uncertainty',
      params: {
        domain: [0, 3],
        list: { kind: 'set', colors: [0x0053d6, 0x65cbf3, 0xffdb13, 0xff7d45] }
      }
    };
  }
  if (m === 'am') {
    // AlphaMissense coloring - B-factors are AM scores in 0-1 range
    // grey (benign) -> orange -> red (pathogenic)
    // NOTE: Molstar's uncertainty theme has `reverse: true`, so colors are reversed
    // Colors listed from high→low B-factor: red, orange, grey, light grey
    return {
      name: 'uncertainty',
      params: {
        domain: [0, 1],
        list: { kind: 'set', colors: [0xd62728, 0xff7d45, 0xbbbbbb, 0xdddddd] }
      }
    };
  }
  if (m === 'dam') {
    // Delta AM coloring - B-factors are delta values in 0-1 range
    // NOTE: Molstar's uncertainty theme has `reverse: true`, so colors are reversed
    return {
      name: 'uncertainty',
      params: {
        domain: [0, 1],
        list: { kind: 'set', colors: [0xd62728, 0xff7d45, 0xbbbbbb, 0xdddddd] }
      }
    };
  }
  if (m === 'aligned') {
    // Aligned vs unaligned: B-factor 1.0 = aligned (green), 0.0 = unaligned (grey)
    // NOTE: Molstar's uncertainty theme has `reverse: true`, so colors are reversed
    // High B-factor (1.0, aligned) → first color, Low (0.0, unaligned) → last color
    return {
      name: 'uncertainty',
      params: {
        domain: [0, 1],
        list: { kind: 'set', colors: [0x43a047, 0xcccccc] }
      }
    };
  }
  if (m === 'domains') {
    // Domain coloring: B-factor 0.0 = no domain (grey), 1.0 = domain (green)
    // NOTE: Molstar's uncertainty theme has `reverse: true`, so colors are reversed
    // High B-factor (1.0, domain) → first color, Low (0.0, no domain) → last color
    return {
      name: 'uncertainty',
      params: {
        domain: [0, 1],
        list: { kind: 'set', colors: [0x2ca02c, 0xcccccc] }
      }
    };
  }
  return {
    name: 'uniform',
    params: { value: 0xcccccc }
  };
}

async function colorBy(mode){
  currentColorMode = mode;
  console.log('Applying color mode:', mode);

  try {
    // Modes that require loading a specific PDB variant with modified B-factors
    const modesPdbMap = {
      'am': () => window.PDB64_AM_BY_MODE && amMode ? window.PDB64_AM_BY_MODE[amMode] : null,
      'plddt': () => window.PDB64_PLDDT || PDB64_FULL,
      'aligned': () => window.PDB64_ALIGNED,
      'domains': () => window.PDB64_DOMAINS,
    };

    if (mode in modesPdbMap) {
      const pdb = modesPdbMap[mode]();
      console.log(`${mode} mode requested, PDB available:`, !!pdb);

      if (pdb) {
        // Debug: check B-factors in the PDB being loaded
        try {
          const pdbText = atob(pdb);
          const lines = pdbText.split('\n').filter(l => l.startsWith('ATOM'));
          const bfactors = lines.slice(0, 20).map(l => parseFloat(l.substring(60, 66).trim()));
          console.log(`${mode} PDB first 20 B-factors:`, bfactors);
          if (mode === 'plddt') {
            // Count bins for pLDDT (bin indices 0-3)
            const allBf = lines.map(l => parseFloat(l.substring(60, 66).trim()));
            const bins = {0: 0, 1: 0, 2: 0, 3: 0};
            allBf.forEach(bf => {
              const bin = Math.round(bf);
              if (bin >= 0 && bin <= 3) bins[bin]++;
            });
            console.log('pLDDT bin distribution:', bins, '(0=orange ≤50, 1=yellow 51-70, 2=cyan 71-90, 3=blue >90)');
          }
        } catch(e) { console.warn('Could not parse PDB for debug:', e); }

        console.log(`Loading ${mode}-colored PDB`);
        await reloadViewerWith(pdb);

        // Wait for structure to be fully loaded
        await new Promise(resolve => setTimeout(resolve, 500));

        // Apply theme after reload
        const theme = themeForColorMode(mode);
        console.log(`Applying ${mode} theme:`, theme);
        await applyColorTheme(theme);

        // Restore selections
        await renderSelections();
        return;
      } else {
        console.warn(`No PDB variant available for mode: ${mode}, falling back to uniform`);
        // Reset to uniform since the requested mode isn't available
        currentColorMode = 'uniform';
        const select = document.getElementById('colorBy');
        if (select) select.value = 'uniform';
      }
    }

    // For uniform mode or if specific PDB not available, just apply theme
    if (!plugin || !structureReady) {
      console.log('Viewer not ready for coloring');
      return;
    }

    const theme = themeForColorMode(mode);
    console.log('Applying theme:', theme.name);
    await applyColorTheme(theme);
    console.log('Color theme applied:', mode);
  } catch(e) {
    console.error('Error applying color theme:', e);
  }
}

async function reloadViewerWith(b64, preserveCamera = false){
  // Save camera state if we want to preserve it
  let cameraSnapshot = null;
  if (preserveCamera && plugin && plugin.canvas3d?.camera) {
    try {
      cameraSnapshot = plugin.canvas3d.camera.getSnapshot();
    } catch(e) { console.warn('Could not save camera state:', e); }
  }

  // If viewer already exists, just load new structure (smoother, no flash)
  if (viewer && plugin) {
    try {
      await plugin.clear();
    } catch(e) {}

    const bytes = Uint8Array.from(atob(b64), c => c.charCodeAt(0));
    const blob = new Blob([bytes], {type:"chemical/x-pdb"});
    const url = URL.createObjectURL(blob);

    try {
      await viewer.loadStructureFromUrl(url, 'pdb');
      structureReady = true;
    } catch(e) {
      console.error('Failed to load structure:', e);
      structureReady = false;
    }

    URL.revokeObjectURL(url);

    // Restore camera state if we saved it
    if (cameraSnapshot && plugin.canvas3d?.camera) {
      try {
        await new Promise(resolve => setTimeout(resolve, 30));
        plugin.canvas3d.camera.setState(cameraSnapshot, 0);
      } catch(e) { console.warn('Could not restore camera state:', e); }
    }
  } else {
    // First load - initialize Molstar
    const host = document.getElementById('viewer');
    host.innerHTML = '';
    viewer = null; plugin = null; structureReady = false;
    await initMolstar();
    await loadPDBfromBase64(b64);
  }
}

async function main(){
  console.log('Initializing report viewer...');
  console.log('PDBe complexes count:', PDBe_COMPLEXES.length);
  initSummarySection();
  
  (DATA.domainsA||[]).forEach(d => { if (d.uid) domByUidA[d.uid] = d; });
  (DATA.domainsB||[]).forEach(d => { if (d.uid) domByUidB[d.uid] = d; });

  document.getElementById('tmScore').textContent = (DATA.tm!=null ? DATA.tm.toFixed(3) : '–');
  document.getElementById('contextTitle').textContent = `Full: ${DATA.g1} × ${DATA.g2}`;
  await reloadViewerWith(PDB64_FULL);

  buildSeq();
  fillDomainTables();
  await fillDomPairs();

  setupPdbeCollapse();
  setupPdbeControls();
  setupAllCollapsibleSections();
  updateSectionVisibility();

  document.getElementById('colorBy').addEventListener('change', (e)=>colorBy(e.target.value), {passive:true});
  document.getElementById('center').addEventListener('click', ()=>{ if(structureReady){ plugin.canvas3d?.requestCameraReset(); }}, {passive:true});
  document.getElementById('lockViewer').addEventListener('click', ()=>{ toggleViewerLock(); }, {passive:true});

  // Set initial lock state (delayed to ensure DOM is ready)
  setTimeout(() => {
    setViewerLocked(true);
  }, 100);
  document.getElementById('backFull').addEventListener('click', async ()=>{
    selection.clear();
    pendingHighlightLoci = null;
    await reloadViewerWith(PDB64_FULL);
    document.getElementById('tmScore').textContent = (DATA.tm!=null ? DATA.tm.toFixed(3) : '–');
    document.getElementById('contextTitle').textContent = `Full: ${DATA.g1} × ${DATA.g2}`;
    Object.values(trackRefs).forEach(track => {
      if (track && track._originalData) track.data = [...track._originalData];
    });
    await initializeHighlightColors();
    setupHoverInterception();
    await renderSelections();
  }, {passive:true});

  // Chain visibility toggles
  const showChainA = document.getElementById('showChainA');
  const showChainB = document.getElementById('showChainB');
  const chainALabel = document.getElementById('chainALabel');
  const chainBLabel = document.getElementById('chainBLabel');

  if (chainALabel) chainALabel.textContent = DATA.g1;
  if (chainBLabel) chainBLabel.textContent = DATA.g2;

  // Track if we're currently processing a chain toggle to prevent double-triggers
  let chainToggleInProgress = false;

  if (showChainA) {
    showChainA.addEventListener('change', async (e) => {
      if (chainToggleInProgress) return;
      chainToggleInProgress = true;

      chainVisible.A = e.target.checked;
      console.log('Chain A toggled:', chainVisible.A, 'Chain B state:', chainVisible.B);

      // Only reset if trying to hide both chains
      if (!chainVisible.A && !chainVisible.B) {
        console.log('Both hidden - resetting to show both');
        chainVisible.A = true;
        chainVisible.B = true;
        showChainA.checked = true;
        if (showChainB) showChainB.checked = true;
      }

      await applyChainVisibility();
      chainToggleInProgress = false;
    }, {passive: true});
  }

  if (showChainB) {
    showChainB.addEventListener('change', async (e) => {
      if (chainToggleInProgress) return;
      chainToggleInProgress = true;

      chainVisible.B = e.target.checked;
      console.log('Chain B toggled:', chainVisible.B, 'Chain A state:', chainVisible.A);

      // Only reset if trying to hide both chains
      if (!chainVisible.A && !chainVisible.B) {
        console.log('Both hidden - resetting to show both');
        chainVisible.A = true;
        chainVisible.B = true;
        if (showChainA) showChainA.checked = true;
        showChainB.checked = true;
      }

      await applyChainVisibility();
      chainToggleInProgress = false;
    }, {passive: true});
  }

  // Druggability filter
  const drugFilter = document.getElementById('druggabilityFilter');
  if (drugFilter) {
    drugFilter.value = 'medium+'; // Default
    drugFilter.addEventListener('change', (e) => {
      druggabilityFilter = e.target.value;
      fillDomainTables(); // Refresh tables with new filter
      applyCavityFilter(); // Also filter Nightingale cavity tracks
    }, {passive: true});
  }

  const amModeSel = document.getElementById('amMode');
  if (amModeSel && Array.isArray(AM_MODES) && AM_MODES.length) {
    amModeSel.innerHTML = '';
    AM_MODES.forEach(m => {
      const opt = document.createElement('option');
      opt.value = m;
      if (m === 'raw') opt.textContent = 'Raw';
      else if (m === 'percentile') opt.textContent = 'Percentile (per protein)';
      else if (m === 'minmax') opt.textContent = 'Min-max (per protein)';
      else if (m === 'zscore') opt.textContent = 'Z-score (logistic)';
      else opt.textContent = m;
      amModeSel.appendChild(opt);
    });
    amModeSel.value = amMode;
    amModeSel.addEventListener('change', (e)=>{
      applyAmMode(e.target.value);
    }, {passive:true});
  }

  // Alignment method toggle
  const alignMethodSel = document.getElementById('alignmentMethod');
  if (alignMethodSel) {
    alignMethodSel.addEventListener('change', async (e) => {
      const method = e.target.value;
      await switchAlignmentMethod(method);
    }, {passive: true});
  }

  setTimeout(async () => {
    await initializeHighlightColors();
    setupHoverInterception();
  }, 1500);
  
  console.log('Report viewer initialized');
}

