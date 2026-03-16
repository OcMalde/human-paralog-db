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

// Get pair ID from URL (let: reassigned in inline mode)
const urlParams = new URLSearchParams(window.location.search);
let PAIR_ID = urlParams.get('pair');

// Family index cache (loaded once)
let FAMILY_INDEX = null;

// These will be populated from API
let DATA = null;
let SUMMARY = null;
let PDB64_FULL = "";
let PLMA_DATA = null;

// PDBe highlight state
let isProteinHighlighted = false;

// Set true when family has only 2 genes — prevents updateSectionVisibility from re-showing PLMA
let isFamily2 = false;

// Lazy section init via IntersectionObserver
function registerLazySection(elementId, initFn) {
  const el = document.getElementById(elementId);
  if (!el || !('IntersectionObserver' in window)) { initFn(); return; }
  const obs = new IntersectionObserver((entries) => {
    if (entries[0].isIntersecting) {
      obs.disconnect();
      initFn();
    }
  }, { rootMargin: '300px' });
  obs.observe(el);
}

// Block Molstar volume server
(function() {
    const block = ['molstarvolseg.ncbr.muni.cz', 'localhost:9000'];
    const _f = window.fetch.bind(window);
    window.fetch = (u, i) => {
        try { if (block.some(h => (typeof u === 'string' ? u : u?.url || '').includes(h))) return Promise.resolve(new Response('{"items":[]}', {status:200})); } catch(e){}
        return _f(u, i);
    };
})();

// Inline mode: load data from window.__INLINE__ (self-contained HTML reports)
async function loadInlineData() {
    const I = window.__INLINE__;
    DATA = I.DATA;
    SUMMARY = I.SUMMARY || { gene1: {}, gene2: {}, pair: {}, conservation: {}, boxplots: {} };
    PLMA_DATA = I.PLMA || null;
    PAIR_ID = DATA.PAIR;

    // PDB variants
    const variants = I.PDB || {};
    PDB64_FULL = variants.pdb64_full || "";
    window.PDB64_A = variants.pdb64_a || PDB64_FULL;
    window.PDB64_B = variants.pdb64_b || PDB64_FULL;
    window.PDB64_AM_BY_MODE = {};
    for (const key in variants) {
        if (key.startsWith('pdb64_am_')) {
            window.PDB64_AM_BY_MODE[key.replace('pdb64_am_', '')] = variants[key];
        }
    }
    window.PDB64_PLDDT = variants.pdb64_plddt || PDB64_FULL;
    window.PDB64_ALIGNED = variants.pdb64_aligned || null;
    window.PDB64_DOMAINS = variants.pdb64_domains || null;

    // DATA-dependent variables
    AM_MODES = DATA.amModes || ['raw'];
    PDBe_COMPLEXES = DATA.pdbeComplexes || [];
    UNIPROT_A = DATA.a1 || '';
    UNIPROT_B = DATA.a2 || '';

    document.title = `${DATA.g1} vs ${DATA.g2}`;
    document.getElementById('titleMain').textContent = `${DATA.g1} ↔ ${DATA.g2}`;
    document.getElementById('titleSub').textContent = `Paralog pair ${DATA.PAIR}`;

    // Family data (inline)
    FULL_FAMILIES = I.FULL_FAMILIES || null;
    FAMILY_INDEX = I.FAMILY_INDEX || null;
    window.__INLINE_INDEX__ = I.INDEX || [];

    await loadFamilyData();
    await main();
    document.getElementById('loadingOverlay').style.display = 'none';
}

// Main initialization - loads data then runs notebook code
async function loadDataAndInit() {
    // Support inline mode (self-contained HTML reports)
    if (window.__INLINE__) return loadInlineData();

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
    let allPairsIndex = [];
    if (window.__INLINE_INDEX__) {
      allPairsIndex = window.__INLINE_INDEX__;
    } else {
      const indexResp = await fetch(`${DATA_BASE}/index.json`);
      allPairsIndex = indexResp.ok ? await indexResp.json() : [];
    }
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

    // 2-gene families: hide both family sections, show a simple note
    if (familyGenes.size <= 2) {
      // familyNav stays hidden (display:none by default — do NOT show it)
      isFamily2 = true; // prevent updateSectionVisibility from re-showing PLMA
      // Hide PLMA section too
      const plmaSection = document.getElementById('familyFeaturesSection');
      if (plmaSection) plmaSection.classList.add('section-hidden');
      // Show a simple message in the family group body
      const groupBody = document.getElementById('familyGroupBody');
      if (groupBody) {
        const note = document.createElement('p');
        note.style.cssText = 'padding:14px 18px;color:#666;font-size:13px;background:#f9f9f9;border-radius:6px;margin:12px 0;border:1px solid #e8e8e8;';
        note.innerHTML = `This pair forms a <strong>family of 2</strong> — phylogenetic tree and family alignment are not available.`;
        groupBody.insertBefore(note, groupBody.firstChild);
      }
      // Grey out sidebar nav links for tree and PLMA
      document.querySelectorAll('.sidebar-nav a[href="#familyNav"], .sidebar-nav a[href="#familyFeaturesSection"]')
        .forEach(a => { a.style.opacity = '0.4'; a.style.pointerEvents = 'none'; a.title = 'Family of 2 — not applicable'; });
      return;
    }

    // Show family sections (only reached when family size > 2)
    familyNav.style.display = 'block';
    const familySubtitle = document.getElementById('familySubtitle');
    if (familySubtitle) {
      familySubtitle.textContent = `${familyGenes.size} genes in family · ${pairsWithReports} pairs with reports`;
    }

    // Initialize constellation immediately, defer tree to not block page load
    setupFamilyViewToggle();
    initFamilyConstellation();
    requestAnimationFrame(() => renderFamilyTree());

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
  // Clamp to the actual container width (avoids overflow when section just expanded)
  const availableWidth = Math.min(container.offsetWidth, container.clientWidth, window.innerWidth - 100);
  const displayWidth = Math.min(Math.max(availableWidth - 20, 300), 1200);
  const displayHeight = Math.min(560, Math.max(400, displayWidth * 0.5));

  // Set canvas size accounting for device pixel ratio (minimum 2x for crisp text)
  canvas.width = displayWidth * dpr;
  canvas.height = displayHeight * dpr;
  canvas.style.width = '100%';
  canvas.style.maxWidth = displayWidth + 'px';
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

// ============= Family Phylogenetic Tree (UPGMA from identity matrix) =============

// Tree subtree helpers
let treeSubtreeOnly = true;

function countTreeLeaves(node) {
  return node.children.length === 0 ? 1 : node.children.reduce((s, c) => s + countTreeLeaves(c), 0);
}

function subtreeContains(node, geneName) {
  if (node.children.length === 0) return node.name === geneName;
  return node.children.some(c => subtreeContains(c, geneName));
}

function findLCA(node, g1, g2) {
  if (!subtreeContains(node, g1) || !subtreeContains(node, g2)) return null;
  for (const child of node.children) {
    const lca = findLCA(child, g1, g2);
    if (lca) return lca;
  }
  return node;
}

function buildUPGMATree(genes, identities) {
  // Build a UPGMA tree from a pairwise identity matrix.
  // genes: array of gene names
  // identities: { gene: { otherGene: identity(0-100) } }
  // Returns a tree: { name, branchLength, children[], size }
  if (genes.length === 0) return null;
  if (genes.length === 1) return { name: genes[0], branchLength: 0, children: [], size: 1 };

  // Build distance matrix (distance = 100 - identity)
  const n = genes.length;
  const dist = Array.from({ length: n }, () => new Float64Array(n));
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const id = (identities[genes[i]] && identities[genes[i]][genes[j]]) ||
                 (identities[genes[j]] && identities[genes[j]][genes[i]]) || 0;
      const d = 100 - id;
      dist[i][j] = d;
      dist[j][i] = d;
    }
  }

  // Clusters: each starts as a leaf
  let clusters = genes.map((g, i) => ({
    node: { name: g, branchLength: 0, children: [], size: 1 },
    indices: [i],
    height: 0,
  }));

  while (clusters.length > 2) {
    // Find closest pair
    let minD = Infinity, mi = -1, mj = -1;
    for (let i = 0; i < clusters.length; i++) {
      for (let j = i + 1; j < clusters.length; j++) {
        let sum = 0, count = 0;
        for (const a of clusters[i].indices) {
          for (const b of clusters[j].indices) {
            sum += dist[a][b];
            count++;
          }
        }
        const avgD = sum / count;
        if (avgD < minD) { minD = avgD; mi = i; mj = j; }
      }
    }

    const newHeight = minD / 2;
    const ci = clusters[mi], cj = clusters[mj];
    ci.node.branchLength = newHeight - ci.height;
    cj.node.branchLength = newHeight - cj.height;

    const merged = {
      node: { name: '', branchLength: 0, children: [ci.node, cj.node], size: ci.node.size + cj.node.size },
      indices: ci.indices.concat(cj.indices),
      height: newHeight,
    };

    // Remove j first (higher index), then i
    clusters.splice(mj, 1);
    clusters.splice(mi, 1);
    clusters.push(merged);
  }

  // Final merge of last two clusters
  if (clusters.length === 2) {
    let sum = 0, count = 0;
    for (const a of clusters[0].indices) {
      for (const b of clusters[1].indices) {
        sum += dist[a][b];
        count++;
      }
    }
    const finalH = (sum / count) / 2;
    clusters[0].node.branchLength = finalH - clusters[0].height;
    clusters[1].node.branchLength = finalH - clusters[1].height;
    return { name: '', branchLength: 0, children: [clusters[0].node, clusters[1].node], size: clusters[0].node.size + clusters[1].node.size };
  }
  return clusters[0].node;
}

function layoutTree(root) {
  // Fixed-angle slanted cladogram: every edge at ±30° from horizontal,
  // every fork exactly 60°. Spine (last child) continues the parent
  // diagonal unbroken; non-spine children branch at 60°.
  // All leaves right-aligned. Positions computed bottom-up from leaves.

  const HALF = Math.PI / 6; // 30°
  const cosA = Math.cos(HALF); // ≈ 0.866
  const sinA = Math.sin(HALF); // 0.5

  function countLeaves(n) {
    return n.children.length === 0 ? 1
      : n.children.reduce((s, c) => s + countLeaves(c), 0);
  }
  const totalLeaves = countLeaves(root);
  if (totalLeaves <= 1) {
    root.x = 0; root.y = 0;
    return { totalLeaves, treeWidth: 0 };
  }

  // Sort children at each node: most leaves last (= spine) for balance
  function sortChildren(n) {
    if (n.children.length > 0) {
      n.children.forEach(sortChildren);
      n.children.sort((a, b) => countLeaves(a) - countLeaves(b));
    }
  }
  sortChildren(root);

  // Binarize polytomies: cascade extra non-spine children
  function binarize(n) {
    n.children.forEach(binarize);
    while (n.children.length > 2) {
      const spine = n.children.pop();           // largest (spine)
      const last = n.children.pop();            // 2nd from end
      const secondLast = n.children.pop();      // 3rd from end
      const virt = { name: '', branchLength: 0, children: [secondLast, last], _virt: true };
      n.children.push(virt);
      n.children.push(spine);
    }
  }
  binarize(root);

  // Assign geometric leaf order (top → bottom).
  // dir=+1: spine goes down-right, non-spine goes up-right
  // dir=-1: spine goes up-right, non-spine goes down-right
  let leafIdx = 0;
  function assignOrder(n, dir) {
    if (n.children.length === 0) { n._li = leafIdx++; return; }
    const spine = n.children[n.children.length - 1];
    const other = n.children[0];
    if (dir === 1) { assignOrder(other, -1); assignOrder(spine, 1); }
    else           { assignOrder(spine, -1); assignOrder(other, 1); }
  }
  assignOrder(root, 1);

  // Set leaf pixel positions (unit space: y = leafIndex, x = 0)
  function initLeaves(n) {
    if (n.children.length === 0) { n.x = 0; n.y = n._li; return; }
    n.children.forEach(initLeaves);
  }
  initLeaves(root);

  // Bottom-up: compute internal node positions from children
  function positionUp(n, dir) {
    if (n.children.length === 0) return;
    const spine = n.children[n.children.length - 1];
    const other = n.children[0];
    positionUp(spine, dir);
    positionUp(other, -dir);
    // Solve for node (nx,ny) such that:
    //   edge to spine at angle dir*30°, edge to other at -dir*30°
    const sumL  = (spine.y - other.y) / (dir * sinA);
    const diffL = (spine.x - other.x) / cosA;
    const Ls = (sumL + diffL) / 2;
    n.y = spine.y - Ls * dir * sinA;
    n.x = spine.x - Ls * cosA;
  }
  positionUp(root, 1);

  // Shift so root.x = 0
  const rx = root.x;
  function shiftX(n) {
    n.x -= rx;
    if (n.children.length > 0) n.children.forEach(shiftX);
  }
  shiftX(root);

  // Compute tree width (max leaf x)
  let maxX = 0;
  function findMaxX(n) {
    if (n.x > maxX) maxX = n.x;
    if (n.children.length > 0) n.children.forEach(findMaxX);
  }
  findMaxX(root);

  return { totalLeaves, treeWidth: maxX };
}

function parseNewick(nwk) {
  let i = 0;
  function parse() {
    const node = { name: '', branchLength: 0, children: [] };
    if (i < nwk.length && nwk[i] === '(') {
      i++;
      while (true) {
        node.children.push(parse());
        if (i >= nwk.length || nwk[i] !== ',') break;
        i++;
      }
      if (i < nwk.length && nwk[i] === ')') i++;
    }
    let name = '';
    while (i < nwk.length && !':,);'.includes(nwk[i])) name += nwk[i++];
    node.name = name.trim();
    if (i < nwk.length && nwk[i] === ':') {
      i++;
      let bl = '';
      while (i < nwk.length && !',);'.includes(nwk[i])) bl += nwk[i++];
      node.branchLength = parseFloat(bl) || 0;
    }
    return node;
  }
  return parse();
}

function renderFamilyTree() {
  const svg = document.getElementById('familyTreeSvg');
  const sourceEl = document.getElementById('treeSourceText');
  if (!svg) return;

  const g1 = DATA ? DATA.g1 : '', g2 = DATA ? DATA.g2 : '';

  // Decide tree source: Compara trees from DATA, or UPGMA fallback
  let trees = []; // Array of tree roots to render
  let treeSource = 'upgma';

  if (DATA && DATA.comparaTrees && DATA.comparaTrees.length > 0) {
    // Parse Compara Newick trees
    for (const ct of DATA.comparaTrees) {
      if (ct.newick) {
        const parsed = parseNewick(ct.newick);
        if (parsed) trees.push(parsed);
      }
    }
    if (trees.length > 0) {
      treeSource = DATA.comparaSource || 'compara';
    }
  }

  // Fallback to UPGMA if no Compara trees (cap at 40 genes for performance)
  if (trees.length === 0 && constellationState && constellationState.allGenes && constellationState.allGenes.length >= 3) {
    let genes = constellationState.allGenes;
    if (genes.length > 40) {
      // Keep query genes + closest by identity
      const g1i = genes.indexOf(g1), g2i = genes.indexOf(g2);
      const scored = genes.map((g, i) => {
        if (i === g1i || i === g2i) return { g, score: Infinity };
        const gd = constellationState.geneData[g];
        const id1 = gd?.identities?.[g1] || 0, id2 = gd?.identities?.[g2] || 0;
        return { g, score: Math.max(id1, id2) };
      });
      scored.sort((a, b) => b.score - a.score);
      genes = scored.slice(0, 40).map(s => s.g);
    }
    const identities = {};
    for (const gene of genes) {
      const gd = constellationState.geneData[gene];
      if (gd && gd.identities) identities[gene] = gd.identities;
    }
    const upgma = buildUPGMATree(genes, identities);
    if (upgma) {
      trees = [upgma];
      treeSource = 'upgma';
    }
  }

  if (trees.length === 0) {
    // No tree at all — hide tree view
    const tc = document.getElementById('familyTreeContainer');
    const th = document.getElementById('familyTreeHelp');
    if (tc) tc.style.display = 'none';
    if (th) th.style.display = 'none';
    const toggle = document.getElementById('familyViewToggle');
    if (toggle) toggle.style.display = 'none';
    const nc = document.getElementById('constellationContainer');
    const nh = document.getElementById('constellationHelp');
    if (nc) nc.style.display = 'flex';
    if (nh) nh.style.display = 'block';
    return;
  }

  // Subtree mode: trim to LCA of g1+g2 if it covers <90% of leaves
  const totalLeavesAll = trees.reduce((s, t) => s + countTreeLeaves(t), 0);
  let showingSubtree = false;
  if (treeSubtreeOnly && trees.length >= 1) {
    let lcaSubtree = null;
    for (const tree of trees) {
      const lca = findLCA(tree, g1, g2);
      if (lca && lca !== tree) { lcaSubtree = lca; break; }
    }
    if (lcaSubtree) {
      const lcaLeaves = countTreeLeaves(lcaSubtree);
      if (lcaLeaves < totalLeavesAll * 0.9) {
        const otherCount = totalLeavesAll - lcaLeaves;
        const stub = { name: `(${otherCount} more…)`, branchLength: 1, children: [], _stub: true };
        trees = [{ name: '', branchLength: 0, children: [stub, lcaSubtree] }];
        showingSubtree = true;
      }
    }
  }

  // Update source text
  if (sourceEl) {
    if (treeSource === 'compara') {
      sourceEl.innerHTML = 'Phylogenetic tree from <strong>Ensembl Compara</strong> (human paralogs only).';
    } else if (treeSource === 'compara_split') {
      sourceEl.innerHTML = 'Separate <strong>Ensembl Compara</strong> sub-trees (genes are in different Compara families).';
    } else {
      sourceEl.innerHTML = 'Tree reconstructed by <strong>UPGMA</strong> from pairwise sequence identity.';
    }
  }

  // Layout all trees and compute total dimensions
  const layouts = trees.map(t => layoutTree(t));
  const splitGap = trees.length > 1 ? 30 : 0;
  let totalLeafRows = layouts.reduce((s, l) => s + l.totalLeaves, 0);
  const rh = totalLeafRows <= 8 ? 22 : (totalLeafRows <= 20 ? 16 : 12);

  // Compute label width from longest gene name
  const allNames = [];
  function collectNames(n) { if (n.children.length === 0 && !n._virt && !n._stub) allNames.push(n.name); else n.children.forEach(collectNames); }
  trees.forEach(collectNames);
  const maxLabelLen = Math.max(...allNames.map(n => n.length), 5);
  const labelW = maxLabelLen * 7 + 16;

  // Width from tree geometry (uniform scale preserves 60° angles)
  const maxTreeW = Math.max(...layouts.map(l => l.treeWidth), 1);
  const plotW = maxTreeW * rh;

  const margin = { top: 12, right: labelW, bottom: 10, left: 10 };
  const treeContentH = totalLeafRows * rh + splitGap * (trees.length - 1);
  const contentH = Math.max(80, treeContentH + margin.top + margin.bottom);
  const svgWidth = Math.max(margin.left + plotW + margin.right, 300);
  const svgHeight = contentH;

  let svgContent = '';
  let yOffset = margin.top;

  for (let ti = 0; ti < trees.length; ti++) {
    const tree = trees[ti];
    const { totalLeaves, treeWidth } = layouts[ti];
    const treeH = totalLeaves * rh;

    // Uniform scaling: preserves the 60° fork angles exactly
    const px = (n) => margin.left + n.x * rh;
    const py = (n) => yOffset + n.y * rh + rh / 2;

    // All edges are straight lines at ±30° — no polylines needed
    function drawBranches(node) {
      if (node.children.length === 0) return;
      const x1 = px(node), y1 = py(node);
      node.children.forEach(c => {
        svgContent += `<line class="tree-branch" x1="${x1}" y1="${y1}" x2="${px(c)}" y2="${py(c)}"/>`;
        drawBranches(c);
      });
    }

    function drawNodes(node) {
      if (node._virt) { node.children.forEach(drawNodes); return; }
      const cx = px(node), cy = py(node);
      if (node.children.length === 0) {
        if (node._stub) {
          svgContent += `<circle class="tree-node other" cx="${cx}" cy="${cy}" r="2.5" opacity="0.4"/>`;
          svgContent += `<text class="tree-label" x="${cx + 6}" y="${cy + 3.5}" fill="#999" font-style="italic">${node.name}</text>`;
        } else {
          const isA = node.name === g1, isB = node.name === g2;
          const cls = isA ? 'gene-a' : isB ? 'gene-b' : 'other';
          svgContent += `<circle class="tree-node ${cls}" cx="${cx}" cy="${cy}" r="2.5"/>`;
          svgContent += `<text class="tree-label ${cls}" x="${cx + 6}" y="${cy + 3.5}" data-gene="${node.name}">${node.name}</text>`;
        }
      } else {
        node.children.forEach(drawNodes);
      }
    }

    drawBranches(tree);
    drawNodes(tree);
    yOffset += treeH;

    // Draw separator between split trees
    if (ti < trees.length - 1) {
      const sepY = yOffset + splitGap / 2;
      svgContent += `<line x1="${margin.left}" y1="${sepY}" x2="${svgWidth - margin.right}" y2="${sepY}" stroke="#ccc" stroke-width="1" stroke-dasharray="6,4"/>`;
      svgContent += `<text x="${svgWidth / 2}" y="${sepY - 6}" text-anchor="middle" fill="#aaa" font-size="10" font-style="italic">separate Compara trees</text>`;
      yOffset += splitGap;
    }
  }

  svg.setAttribute('viewBox', `0 0 ${svgWidth} ${svgHeight}`);
  svg.setAttribute('preserveAspectRatio', 'xMidYMin meet');
  svg.removeAttribute('width');
  svg.removeAttribute('height');
  svg.style.width = '100%';
  svg.style.height = 'auto';
  svg.style.maxWidth = '100%';
  svg.innerHTML = svgContent;

  // Wire click on gene labels
  svg.querySelectorAll('.tree-label').forEach(el => {
    el.addEventListener('click', () => {
      const gene = el.dataset.gene;
      if (!gene) return;
      const tryIds = [
        `${g1}_${gene}`, `${gene}_${g1}`,
        `${g2}_${gene}`, `${gene}_${g2}`,
      ];
      for (const pid of tryIds) {
        if (constellationState && constellationState.pairData && constellationState.pairData[pid]) {
          window.location.href = `?pair=${pid}`;
          return;
        }
      }
    });
  });

  // Show/hide "expand tree" toggle button
  let expandBtn = document.getElementById('treeExpandBtn');
  if (totalLeavesAll > 3) {
    if (!expandBtn) {
      expandBtn = document.createElement('button');
      expandBtn.id = 'treeExpandBtn';
      expandBtn.style.cssText = 'margin-top:6px;font-size:11px;padding:3px 10px;border:1px solid #bbb;border-radius:4px;background:#f5f5f5;cursor:pointer;color:#555';
      const help = document.getElementById('familyTreeHelp');
      if (help) help.insertBefore(expandBtn, help.firstChild);
    }
    expandBtn.textContent = showingSubtree ? `Show full tree (${totalLeavesAll} genes)` : 'Show subtree only';
    expandBtn.onclick = () => {
      treeSubtreeOnly = !treeSubtreeOnly;
      renderFamilyTree();
    };
    expandBtn.style.display = '';
  } else if (expandBtn) {
    expandBtn.style.display = 'none';
  }
}

function setupFamilyViewToggle() {
  const btns = document.querySelectorAll('#familyViewToggle .fv-btn');
  const treeContainer = document.getElementById('familyTreeContainer');
  const treeHelp = document.getElementById('familyTreeHelp');
  const netContainer = document.getElementById('constellationContainer');
  const netHelp = document.getElementById('constellationHelp');

  btns.forEach(btn => {
    btn.addEventListener('click', () => {
      btns.forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      const view = btn.dataset.view;
      if (view === 'tree') {
        if (treeContainer) treeContainer.style.display = 'flex';
        if (treeHelp) treeHelp.style.display = 'block';
        if (netContainer) netContainer.style.display = 'none';
        if (netHelp) netHelp.style.display = 'none';
      } else {
        if (treeContainer) treeContainer.style.display = 'none';
        if (treeHelp) treeHelp.style.display = 'none';
        if (netContainer) netContainer.style.display = 'flex';
        if (netHelp) netHelp.style.display = 'block';
        initFamilyConstellation();
      }
    });
  });
}

// ============= End Family Tree =============

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
  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, width, height);

  // Draw orbit circles (identity thresholds) - including 0% and 100%
  const thresholds = [0, 0.2, 0.4, 0.6, 0.8, 1.0];
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
  const searchIndexPromise = window.__INLINE_INDEX__
    ? Promise.resolve(window.__INLINE_INDEX__)
    : fetch(`${DATA_BASE}/index.json`).then(resp => resp.json());
  searchIndexPromise.then(pairs => {
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

// Track group toggle state (persists across buildSeq rebuilds)
const trackGroupState = { default: true, structure: false, conservation: false, druggability: false };
let trackGroupRows = {}; // group name -> [row elements]
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
  // DxD header labels
  const dxdHA = document.getElementById('dxdHA');
  const dxdHB = document.getElementById('dxdHB');
  if (dxdHA) dxdHA.textContent = `${g1} domain`;
  if (dxdHB) dxdHB.textContent = `${g2} domain`;
  
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

  // Initialize new sections (non-canvas ones are immediate)
  initProteinDescriptions();
  initSlFunctionalSection();
  // Canvas-heavy sections are deferred to lazy loading (registered in main())

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

// Convert PubMed IDs in text to clickable links
function linkifyPubMed(text) {
  // Match patterns like "PubMed:12345678" or "(PubMed:12345678)"
  return text.replace(/PubMed:(\d+)/g,
    '<a href="https://pubmed.ncbi.nlm.nih.gov/$1" target="_blank" style="color:#1565c0">PubMed:$1</a>');
}

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

    funcEl.innerHTML = linkifyPubMed(truncated);
    funcEl.classList.add('truncated');
    funcEl.dataset.fullText = fullText;
    funcEl.dataset.truncatedText = truncated;
    toggleBtn.style.display = 'inline-block';
    toggleBtn.textContent = 'Show more';

    toggleBtn.addEventListener('click', () => {
      const isExpanded = funcEl.classList.contains('expanded');
      if (isExpanded) {
        funcEl.innerHTML = linkifyPubMed(funcEl.dataset.truncatedText);
        funcEl.classList.remove('expanded');
        funcEl.classList.add('truncated');
        toggleBtn.textContent = 'Show more';
      } else {
        funcEl.innerHTML = linkifyPubMed(funcEl.dataset.fullText);
        funcEl.classList.remove('truncated');
        funcEl.classList.add('expanded');
        toggleBtn.textContent = 'Show less';
      }
    }, { passive: true });
  } else {
    funcEl.innerHTML = linkifyPubMed(fullText);
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

  // Known drugs — now rendered in unified table section
  initKnownDrugsSection();
}

// Track selected drug globally (by drugId)
let _selectedDrugId = null;

function initKnownDrugsSection() {
  const container = document.getElementById('knownDrugsTable');
  if (!container) return;

  const gene1 = SUMMARY.gene1 || {};
  const gene2 = SUMMARY.gene2 || {};
  const g1 = DATA.g1 || gene1.symbol || 'Gene A';
  const g2 = DATA.g2 || gene2.symbol || 'Gene B';
  const drugsA = gene1.known_drugs || [];
  const drugsB = gene2.known_drugs || [];

  if (drugsA.length === 0 && drugsB.length === 0) {
    container.innerHTML = '<p style="color:#888;font-size:13px;margin:0;">No known drugs for either gene.</p>';
    return;
  }

  // Build unified drug list with target info
  const drugMap = new Map();
  for (const d of drugsA) {
    drugMap.set(d.drugId, { ...d, targetsA: true, targetsB: false });
  }
  for (const d of drugsB) {
    if (drugMap.has(d.drugId)) {
      drugMap.get(d.drugId).targetsB = true;
    } else {
      drugMap.set(d.drugId, { ...d, targetsA: false, targetsB: true });
    }
  }
  const allDrugs = [...drugMap.values()];

  // Sort: shared first, then by phase desc, then name
  allDrugs.sort((a, b) => {
    const sharedA = a.targetsA && a.targetsB ? 1 : 0;
    const sharedB = b.targetsA && b.targetsB ? 1 : 0;
    if (sharedB !== sharedA) return sharedB - sharedA;
    if ((b.phase || 0) !== (a.phase || 0)) return (b.phase || 0) - (a.phase || 0);
    return (a.name || '').localeCompare(b.name || '');
  });

  const MAX_VISIBLE = 8;
  const needsExpand = allDrugs.length > MAX_VISIBLE;

  // Phase pill helper
  function phasePill(phase) {
    const p = phase || 0;
    const label = p > 0 ? `Phase ${p}` : 'Unknown';
    return `<span class="phase-pill phase-${p}">${label}</span>`;
  }

  // Target pill
  function targetPill(d) {
    const isShared = d.targetsA && d.targetsB;
    if (isShared) return `<span class="target-pill target-both">Both</span>`;
    return `<span class="target-pill target-single">${d.targetsA ? g1 : g2}</span>`;
  }

  // Build table
  let html = '<table class="drug-table"><thead><tr>';
  html += '<th>Drug</th><th>Phase</th><th>Type</th><th>Targets</th><th>Mechanism</th>';
  html += '</tr></thead><tbody>';

  allDrugs.forEach((d, i) => {
    const isShared = d.targetsA && d.targetsB;
    const hidden = needsExpand && i >= MAX_VISIBLE ? ' style="display:none" data-drug-extra="1"' : '';
    const rowClass = isShared ? 'drug-shared' : '';
    const name = d.name || d.drugId || '?';
    const moa = d.moa || '';
    const moaShort = moa.length > 50 ? moa.substring(0, 47) + '...' : moa;
    html += `<tr class="${rowClass}" data-drug-id="${d.drugId || ''}"${hidden}>`;
    html += `<td class="drug-name">${name}</td>`;
    html += `<td>${phasePill(d.phase)}</td>`;
    html += `<td style="font-size:11px;color:#666">${d.type || ''}</td>`;
    html += `<td>${targetPill(d)}</td>`;
    html += `<td style="font-size:11px;color:#777" title="${moa}">${moaShort}</td>`;
    html += '</tr>';
  });
  html += '</tbody></table>';

  if (needsExpand) {
    html += `<button class="drug-expand-btn" id="drugExpandBtn">Show ${allDrugs.length - MAX_VISIBLE} more drugs</button>`;
  }

  // Summary line
  const nShared = allDrugs.filter(d => d.targetsA && d.targetsB).length;
  const nTotal = allDrugs.length;
  html = `<div style="font-size:12px;color:#666;margin-bottom:8px;">${nTotal} drugs total (${drugsA.length} for ${g1}, ${drugsB.length} for ${g2}${nShared > 0 ? `, <strong style="color:#283593">${nShared} shared</strong>` : ''})</div>` + html;

  container.innerHTML = html;

  // Expand button
  if (needsExpand) {
    const btn = document.getElementById('drugExpandBtn');
    if (btn) {
      let expanded = false;
      btn.addEventListener('click', () => {
        expanded = !expanded;
        container.querySelectorAll('[data-drug-extra]').forEach(tr => {
          tr.style.display = expanded ? '' : 'none';
        });
        btn.textContent = expanded ? 'Show fewer' : `Show ${allDrugs.length - MAX_VISIBLE} more drugs`;
      });
    }
  }

  // Click handler for drug rows — highlight matching drugs
  container.querySelectorAll('.drug-name').forEach(cell => {
    cell.addEventListener('click', () => {
      const tr = cell.closest('tr');
      const drugId = tr?.dataset.drugId;
      if (drugId) toggleKnownDrug(drugId);
    });
  });
}

function toggleKnownDrug(drugId) {
  if (!drugId) return;
  _selectedDrugId = (_selectedDrugId === drugId) ? null : drugId;
  // Highlight in drug table
  document.querySelectorAll('.drug-table tr[data-drug-id]').forEach(tr => {
    if (tr.dataset.drugId === _selectedDrugId) {
      tr.style.background = '#bbdefb';
      tr.style.outline = '2px solid #1976d2';
      tr.style.outlineOffset = '-1px';
    } else {
      tr.style.background = '';
      tr.style.outline = '';
    }
  });
}

// Boxplot chart instances for new sections
// (boxplot chart removed - replaced by PLMA alignment)

// Current similarity search mode (struct or seq) and view (0=overview, 1=scale)
let simSearchMode = 'struct';
let simSearchView = 0;
let simSearchQueryGene = 'A'; // 'A' or 'B'
// Store hit-test regions for hover tooltips
let simSearchHitRegions = [];
let simSearchBarHitRegions = [];

function initSimilaritySearchSection() {
  const canvas = document.getElementById('simSearchRankCanvas');
  const barCanvas = document.getElementById('simSearchBarCanvas');
  const modeSelect = document.getElementById('simSearchModeSelect');
  if (!canvas) return;

  // Query gene toggle: show if reverse data exists
  const queryToggle = document.getElementById('simSearchQueryToggle');
  if (queryToggle) {
    const simSearch = SUMMARY?.similarity_search || {};
    const hasReverse = !!(simSearch['rank_struct_reverse'] || simSearch['rank_seq_reverse']);
    if (hasReverse) {
      queryToggle.style.display = '';
      // Set gene labels
      const la = document.getElementById('simQueryLabelA');
      const lb = document.getElementById('simQueryLabelB');
      if (la) la.textContent = SUMMARY?.gene1?.symbol || 'Gene A';
      if (lb) lb.textContent = SUMMARY?.gene2?.symbol || 'Gene B';
      // Wire radios
      queryToggle.querySelectorAll('input[name="simQuery"]').forEach(r => {
        r.addEventListener('change', () => {
          simSearchQueryGene = r.value;
          drawSimSearchRankViz(simSearchMode);
          drawSimSearchBarViz(simSearchMode);
        });
      });
    }
  }

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
    if (btn0) btn0.style.background = idx === 0 ? '#e0f2f1' : '';
    if (btn0) btn0.style.color = idx === 0 ? '#00695c' : '';
    if (btn0) btn0.style.fontWeight = idx === 0 ? '600' : '';
    if (btn1) btn1.style.background = idx === 1 ? '#e0f2f1' : '';
    if (btn1) btn1.style.color = idx === 1 ? '#00695c' : '';
    if (btn1) btn1.style.fontWeight = idx === 1 ? '600' : '';
  }
  if (btn0) btn0.addEventListener('click', () => setView(0));
  if (btn1) btn1.addEventListener('click', () => setView(1));
  setView(1);

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
  const rev = simSearchQueryGene === 'B' ? '_reverse' : '';
  const fwdRank = simSearch['rank' + suffix];
  const revRank = simSearch['rank' + suffix + '_reverse'];
  const fwdSP   = simSearch['selfSP' + suffix];
  const revSP   = simSearch['selfSP' + suffix + '_reverse'];
  const fwdTax  = simSearch['taxid' + suffix];
  const revTax  = simSearch['taxid' + suffix + '_reverse'];
  // Global max = worst across both directions (fixed scale regardless of query direction)
  const globalRankMax   = Math.max(fwdRank?.max_value || 1, revRank?.max_value || 1);
  const globalSPMax     = Math.max(fwdSP?.max_value   || 1, revSP?.max_value   || 1);
  const globalTaxidMax  = Math.max(fwdTax?.max_value  || 1, revTax?.max_value  || 1);
  const curRank = rev ? revRank : fwdRank;
  const curSP   = rev ? revSP   : fwdSP;
  const curTax  = rev ? revTax  : fwdTax;
  return {
    suffix,
    dbFullName: mode === 'struct' ? 'AlphaFold DB' : 'UniProt',
    rankInfo:  curRank  || fwdRank,
    selfSPInfo: curSP   || fwdSP,
    taxidInfo:  curTax  || fwdTax,
    globalRankMax, globalSPMax, globalTaxidMax,
    gene1: simSearchQueryGene === 'B' ? (SUMMARY?.gene2?.symbol || 'B') : (SUMMARY?.gene1?.symbol || 'A'),
    gene2: simSearchQueryGene === 'B' ? (SUMMARY?.gene1?.symbol || 'A') : (SUMMARY?.gene2?.symbol || 'B'),
  };
}

// Compute bar fill % for rank: linear scale capped at p95
function rankLinearPct(value, p95) {
  if (value == null || p95 == null || p95 <= 0) return 0;
  return Math.min(100, (value / p95) * 100);
}

// ====== VIEW 1: Overview (metrics box) ======
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

  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, displayWidth, displayHeight);

  const m = ssGetMetrics(mode);
  const rank = m.rankInfo?.value ?? null;
  const selfSP = m.selfSPInfo?.value ?? null;
  const taxid = m.taxidInfo?.value ?? null;

  const rankP95 = m.rankInfo?.p95_value ?? 500;
  const rankFillPct = rankLinearPct(rank, rankP95);
  const selfSPPct = m.selfSPInfo?.percentile ?? 50;
  const selfSPPctRel = m.selfSPInfo?.percentile_rank_relative ?? selfSPPct;
  const taxidPct = m.taxidInfo?.percentile ?? 50;
  const taxidPctRel = m.taxidInfo?.percentile_rank_relative ?? taxidPct;
  const totalPairs = m.rankInfo?.total_pairs ?? 105107;

  const dbBoxColor = '#f5f1e6';
  const dbBoxStroke = '#e0d7c2';
  const barBgColor = '#eae6dc';
  const barFillColor = '#5f4d2f';
  const barFillRelColor = '#8b7a5e';
  const barBorderColor = '#333';

  const dbCenterX = displayWidth / 2;
  const dbBoxWidth = 320;
  const dbBoxHeight = 220;
  const dbBoxLeft = dbCenterX - dbBoxWidth / 2;
  const dbBoxTop = (displayHeight - dbBoxHeight) / 2;
  const dbBoxRadius = 12;

  // Direction label above box
  ctx.font = '13px -apple-system, sans-serif';
  ctx.fillStyle = '#555';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'bottom';
  ctx.fillText(`Query: ${m.gene1}  →  ${m.gene2}`, dbCenterX, dbBoxTop - 8);

  function drawBar(x, y, width, height, fillPct, fillColor, r) {
    r = r || 5;
    ssRoundRect(ctx, x, y, width, height, r);
    ctx.fillStyle = barBgColor;
    ctx.fill();
    ctx.strokeStyle = barBorderColor;
    ctx.lineWidth = 1;
    ctx.stroke();
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

  function drawSplitBar(x, y, width, pctTop, pctBottom) {
    const bh = 12;
    const gap = 3;
    drawBar(x, y, width, bh, pctTop, barFillColor, 4);
    drawBar(x, y + bh + gap, width, bh, pctBottom, barFillRelColor, 4);
  }

  // Database box
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
  const metricSpacing = 58;
  const labelX = dbBoxLeft + 14;
  const barX = dbBoxLeft + 14;
  const barWidth = dbBoxWidth - 28;

  const metricLabels = [
    'Rank of the paralog',
    'Human proteins ranking better',
    'Species with proteins ranking better',
  ];

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

  // --- selfSP row ---
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

  // --- taxid row ---
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

// ====== VIEW 2: Scale bars (rank 1 ====▼ - - - worst) ======
function drawSimSearchBarViz(mode) {
  const canvas = document.getElementById('simSearchBarCanvas');
  if (!canvas || !SUMMARY) return;

  const ctx = canvas.getContext('2d');
  const dpr = Math.max(2, window.devicePixelRatio || 1);
  simSearchBarHitRegions = [];

  const displayWidth = 680;
  const displayHeight = 310;
  canvas.width = displayWidth * dpr;
  canvas.height = displayHeight * dpr;
  canvas.style.width = displayWidth + 'px';
  canvas.style.height = displayHeight + 'px';
  ctx.scale(dpr, dpr);

  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, displayWidth, displayHeight);

  const m = ssGetMetrics(mode);
  const rank = m.rankInfo?.value ?? null;
  const selfSP = m.selfSPInfo?.value ?? null;
  const taxid = m.taxidInfo?.value ?? null;
  const totalPairs = m.rankInfo?.total_pairs ?? 105107;

  // Use GLOBAL max (worst across both query directions) so scale is stable
  const globalMaxValues = [m.globalRankMax || 1, m.globalSPMax || 1, m.globalTaxidMax || 1];

  const filledColor = '#e8dcc8';
  const filledBorder = '#333';
  const emptyBorder = '#999';
  const markerColor = '#4a6741';  // dark green marker

  const barLabels = [
    'Rank of the paralog',
    'Human proteins ranking better',
    'Species with proteins ranking better',
  ];
  const values = [rank, selfSP, taxid];
  const infos = [m.rankInfo, m.selfSPInfo, m.taxidInfo];

  const padLeft = 16;
  const padRight = 90;  // room for "worst (N)" label with value
  const barTotalWidth = displayWidth - padLeft - padRight;
  const barHeight = 22;
  const barSpacing = 90;
  const startY = 36;

  // Direction label at top
  ctx.font = 'bold 13px -apple-system, sans-serif';
  ctx.fillStyle = '#444';
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.fillText(`Query: ${m.gene1}  →  ${m.gene2}`, padLeft, 6);

  for (let i = 0; i < 3; i++) {
    const val = values[i];
    const globalMax = globalMaxValues[i];
    const info = infos[i];
    const by = startY + i * barSpacing;
    const cy = by + barHeight / 2;

    // Label right-aligned above bar (avoid overlap with triangle)
    ctx.font = '11px -apple-system, sans-serif';
    ctx.fillStyle = '#666';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'bottom';
    ctx.fillText(barLabels[i], padLeft + barTotalWidth, by - 2);

    if (val == null || globalMax <= 0) {
      // Draw empty bar with N/A
      ssRoundRect(ctx, padLeft, by, barTotalWidth, barHeight, 5);
      ctx.fillStyle = '#f0ece2';
      ctx.fill();
      ctx.strokeStyle = '#ccc'; ctx.lineWidth = 1; ctx.stroke();
      continue;
    }

    const fillFrac = Math.min(1, Math.max(0, val / globalMax));
    const markerX = padLeft + fillFrac * barTotalWidth;

    // Full grey background bar
    ssRoundRect(ctx, padLeft, by, barTotalWidth, barHeight, 5);
    ctx.fillStyle = '#f0ece2';
    ctx.fill();
    ctx.strokeStyle = '#ccc';
    ctx.lineWidth = 1;
    ctx.stroke();

    // Filled portion (rank 1 → pair position)
    if (fillFrac > 0) {
      const fw = Math.max(2, fillFrac * barTotalWidth);
      ssRoundRect(ctx, padLeft, by, fw, barHeight, 5);
      ctx.fillStyle = filledColor;
      ctx.fill();
      ctx.strokeStyle = filledBorder;
      ctx.lineWidth = 1;
      ctx.stroke();
    }

    // Empty dashed portion (pair position → worst)
    const emptyStartX = padLeft + Math.max(2, fillFrac * barTotalWidth);
    const emptyW = barTotalWidth - Math.max(2, fillFrac * barTotalWidth);
    if (emptyW > 4) {
      ssRoundRect(ctx, emptyStartX, by, emptyW, barHeight, 5);
      ctx.fillStyle = '#ffffff';
      ctx.fill();
      ctx.setLineDash([4, 4]);
      ctx.strokeStyle = emptyBorder;
      ctx.lineWidth = 1;
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // Marker triangle (pointing down) above the bar at pair's position
    const triH = 10, triW = 10;
    ctx.beginPath();
    ctx.moveTo(markerX, by - 2);
    ctx.lineTo(markerX - triW / 2, by - 2 - triH);
    ctx.lineTo(markerX + triW / 2, by - 2 - triH);
    ctx.closePath();
    ctx.fillStyle = markerColor;
    ctx.fill();

    // Small "1" tick at left edge below bar (axis reference)
    ctx.font = '9px -apple-system, sans-serif';
    ctx.fillStyle = '#aaa';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    ctx.fillText('1', padLeft, by + barHeight + 3);

    // Protein name + rank centered at marker position, below bar
    // Clamp so it doesn't collide with "1" on left or "worst" label on right
    const labelText = `${m.gene2}: ${formatNum(val)}`;
    ctx.font = 'bold 10px -apple-system, sans-serif';
    const labelW = ctx.measureText(labelText).width;
    const markerLabelX = Math.min(
      Math.max(markerX, padLeft + labelW / 2 + 14),
      padLeft + barTotalWidth - labelW / 2
    );
    ctx.fillStyle = markerColor;
    ctx.textAlign = 'center';
    ctx.fillText(labelText, markerLabelX, by + barHeight + 3);

    // "worst (N)" label at right end — always shows global max value
    ctx.font = '9px -apple-system, sans-serif';
    ctx.fillStyle = '#aaa';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    ctx.fillText(`worst (${formatNum(globalMax)})`, padLeft + barTotalWidth + 6, cy);

    // Hit region for tooltip
    const rPos = info?.rank_position ?? null;
    const pctStr = (fillFrac * 100).toFixed(1);
    simSearchBarHitRegions.push({
      x: padLeft, y: by, w: barTotalWidth, h: barHeight,
      tooltip: rPos != null
        ? `${m.gene2}: rank ${formatNum(val)} (${pctStr}% of scale)<br>${formatNum(rPos)} of ${formatNum(totalPairs)} pairs have ≤ ${formatNum(val)}`
        : `${m.gene2}: rank ${formatNum(val)}`,
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

  // "Show all paralogs" toggle — disable if ≤ 2 human members
  const allParalogsCheckbox = document.getElementById('plmaShowAllParalogs');
  if (allParalogsCheckbox) {
    const nHumanTotal = (PLMA_DATA.sequences || []).filter(s => s.is_human !== false).length;
    if (nHumanTotal <= 2) {
      allParalogsCheckbox.disabled = true;
      allParalogsCheckbox.parentElement.style.opacity = '0.4';
      allParalogsCheckbox.parentElement.title = 'Only 2 family members';
    }
    allParalogsCheckbox.addEventListener('change', () => { drawPlmaAlignment(); });
  }

  // Ortholog toggle — grey out if no orthologs present
  const orthoCheckbox = document.getElementById('plmaShowOrthologs');
  if (orthoCheckbox) {
    const hasOrthologs = (PLMA_DATA.sequences || []).some(s => s.is_human === false);
    if (!hasOrthologs) {
      orthoCheckbox.disabled = true;
      orthoCheckbox.parentElement.style.opacity = '0.4';
      orthoCheckbox.parentElement.title = 'No orthologs in this alignment';
    }
    orthoCheckbox.addEventListener('change', () => { drawPlmaAlignment(); });
  }

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

// --- PLMA grouping helpers ---

function plmaJaccard(setA, setB) {
  let inter = 0;
  for (const v of setA) { if (setB.has(v)) inter++; }
  const union = setA.size + setB.size - inter;
  return union > 0 ? inter / union : 0;
}

function buildSimplifiedParalogRows(humanSeqs, blocks, geneASeq, geneBSeq) {
  // Build block-presence sets
  const blockSets = new Map();
  for (const seq of humanSeqs) {
    const myBlocks = new Set();
    for (let bi = 0; bi < blocks.length; bi++) {
      if (blocks[bi].positions[seq.num]) myBlocks.add(bi);
    }
    blockSets.set(seq.num, myBlocks);
  }
  // Pair blocks
  const pairBlocks = new Set();
  for (let bi = 0; bi < blocks.length; bi++) {
    if (blocks[bi].positions[geneASeq] || blocks[bi].positions[geneBSeq]) pairBlocks.add(bi);
  }
  // Filter: keep paralogs sharing ≥ 1 block with pair
  const relevant = humanSeqs.filter(seq => {
    const myBlocks = blockSets.get(seq.num);
    for (const bi of myBlocks) { if (pairBlocks.has(bi)) return true; }
    return false;
  });
  // Group by Jaccard ≥ 0.7 (complete-linkage greedy)
  const THRESHOLD = 0.7;
  const assigned = new Set();
  const groups = [];
  for (const seq of relevant) {
    if (assigned.has(seq.num)) continue;
    const group = [seq];
    assigned.add(seq.num);
    for (const cand of relevant) {
      if (assigned.has(cand.num)) continue;
      let allOk = true;
      const candSet = blockSets.get(cand.num);
      for (const mem of group) {
        if (plmaJaccard(candSet, blockSets.get(mem.num)) < THRESHOLD) { allOk = false; break; }
      }
      if (allOk) { group.push(cand); assigned.add(cand.num); }
    }
    groups.push(group);
  }
  return groups;
}

// --- Main PLMA drawing ---

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

  const showOrthologs = document.getElementById('plmaShowOrthologs')?.checked || false;
  const showAllParalogs = document.getElementById('plmaShowAllParalogs')?.checked || false;

  const seqA = sequences.find(s => s.num === geneASeq);
  const seqB = sequences.find(s => s.num === geneBSeq);
  const humanParalogs = [];
  const orthologs = [];
  for (const s of sequences) {
    if (s.num === geneASeq || s.num === geneBSeq) continue;
    if (s.is_human === false) orthologs.push(s);
    else humanParalogs.push(s);
  }

  // Build displayRows
  const displayRows = [];
  const mkLabel = (s) => s.gene || s.uniprot || `Seq ${s.num}`;

  if (seqA) displayRows.push({
    label: mkLabel(seqA), labelColor: '#92400e',
    isPair: true, isPairA: true, isPairB: false, isGroup: false,
    members: [seqA], seqNums: new Set([seqA.num]),
    font: 'bold 12px -apple-system, sans-serif',
  });
  if (seqB) displayRows.push({
    label: mkLabel(seqB), labelColor: '#5b21b6',
    isPair: true, isPairA: false, isPairB: true, isGroup: false,
    members: [seqB], seqNums: new Set([seqB.num]),
    font: 'bold 12px -apple-system, sans-serif',
  });

  if (showAllParalogs || humanParalogs.length <= 2) {
    for (const s of humanParalogs) {
      displayRows.push({
        label: mkLabel(s), labelColor: '#666',
        isPair: false, isPairA: false, isPairB: false, isGroup: false,
        members: [s], seqNums: new Set([s.num]),
        font: '11px -apple-system, sans-serif',
      });
    }
  } else {
    const groups = buildSimplifiedParalogRows(humanParalogs, blocks, geneASeq, geneBSeq);
    for (const group of groups) {
      const names = group.map(mkLabel);
      const label = names.length === 1 ? names[0] : names.length <= 3 ? names.join(', ') : `${names[0]} +${names.length - 1}`;
      displayRows.push({
        label, labelColor: group.length > 1 ? '#2563eb' : '#666',
        isPair: false, isPairA: false, isPairB: false,
        isGroup: group.length > 1,
        members: group, seqNums: new Set(group.map(s => s.num)),
        font: '11px -apple-system, sans-serif',
      });
    }
  }

  if (showOrthologs) {
    for (const s of orthologs) {
      displayRows.push({
        label: mkLabel(s), labelColor: '#999',
        isPair: false, isPairA: false, isPairB: false, isGroup: false,
        members: [s], seqNums: new Set([s.num]),
        font: '10px -apple-system, sans-serif',
      });
    }
  }

  const nRows = displayRows.length;
  const nBlocks = blocks.length;

  // --- Merge consecutive blocks with same presence+category (simplified view) ---
  const isSimplified = !showAllParalogs && humanParalogs.length > 2;
  function blockFingerprint(bi) {
    const cat = blocks[bi].category;
    const bits = displayRows.map(row =>
      row.isGroup
        ? row.members.map(m => blocks[bi].positions[m.num] ? '1' : '0').join('')
        : (row.members.some(m => blocks[bi].positions[m.num]) ? '1' : '0')
    ).join(',');
    return cat + ':' + bits;
  }
  const mergedGroups = [];
  if (isSimplified && nBlocks > 0) {
    let i = 0;
    while (i < nBlocks) {
      const fp = blockFingerprint(i);
      const grp = [i];
      while (i + 1 < nBlocks && blockFingerprint(i + 1) === fp) { i++; grp.push(i); }
      mergedGroups.push(grp);
      i++;
    }
  } else {
    for (let i = 0; i < nBlocks; i++) mergedGroups.push([i]);
  }
  const nMerged = mergedGroups.length;

  // Dynamic label width
  let labelWidth = 90;
  for (const row of displayRows) {
    ctx.font = row.font;
    const w = ctx.measureText(row.label).width;
    if (w + 16 > labelWidth) labelWidth = Math.min(200, Math.ceil(w + 16));
  }

  const padRight = 24;
  const padTop = 20;
  const trackHeight = nRows <= 6 ? 22 : (nRows <= 15 ? 14 : 10);
  const trackGap = nRows <= 6 ? 10 : (nRows <= 15 ? 6 : 3);
  const pairTrackHeight = trackHeight + 6;

  const scrollWrapper = document.getElementById('plmaScrollWrapper');
  const containerWidth = Math.max(700, scrollWrapper?.clientWidth || canvas.parentElement?.clientWidth || 700);

  const blockMaxLen = blocks.map(b => {
    let mx = 0;
    for (const pos of Object.values(b.positions)) mx = Math.max(mx, pos.length || 0);
    return mx;
  });
  const groupMaxLen = mergedGroups.map(g => g.reduce((s, bi) => s + blockMaxLen[bi], 0));
  const totalBlockAA = groupMaxLen.reduce((a, b) => a + b, 0) || 1;

  const nGaps = Math.max(0, nMerged - 1);
  const gapWidthBase = nMerged <= 10 ? 14 : (nMerged <= 50 ? 8 : 4);
  const minPxPerAA = 0.6;
  const minBlockArea = totalBlockAA * minPxPerAA;
  const minNeededWidth = labelWidth + padRight + minBlockArea + nGaps * gapWidthBase;
  let displayWidth = Math.max(containerWidth, minNeededWidth);
  const trackAreaWidth = displayWidth - labelWidth - padRight;

  const totalGapWidth = nGaps * gapWidthBase;
  const blockAreaWidth = Math.max(trackAreaWidth - totalGapWidth, trackAreaWidth * 0.5);
  const gapWidth = nGaps > 0 ? (trackAreaWidth - blockAreaWidth) / nGaps : 0;
  const minBlockPx = 3;

  const groupColX = [];
  const groupColW = [];
  let curX = labelWidth;
  for (let gi = 0; gi < nMerged; gi++) {
    groupColX.push(curX);
    const w = Math.max(minBlockPx, (groupMaxLen[gi] / totalBlockAA) * blockAreaWidth);
    groupColW.push(w);
    curX += w;
    if (gi < nMerged - 1) curX += gapWidth;
  }
  // Ensure canvas is wide enough for all blocks + right padding + stroke buffer
  const actualNeeded = curX + padRight + 2;
  if (actualNeeded > displayWidth) {
    displayWidth = Math.ceil(actualNeeded);
  }

  // Canvas height
  const totalTrackHeight = pairTrackHeight * 2 + (trackGap * 2) +
    Math.max(0, nRows - 2) * (trackHeight + trackGap) + padTop + 30;
  const displayHeight = Math.max(180, totalTrackHeight);

  canvas.width = displayWidth * dpr;
  canvas.height = displayHeight * dpr;
  canvas.style.width = displayWidth + 'px';
  canvas.style.height = displayHeight + 'px';
  const inner = document.getElementById('plmaCanvasInner');
  if (inner) inner.style.minWidth = displayWidth + 'px';
  ctx.scale(dpr, dpr);

  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, displayWidth, displayHeight);

  const catColors = {
    specific_a:         '#EF5350',
    specific_b:         '#AB47BC',
    pair_exclusive:     '#26A69A',
    a_with_family:      '#FF7043',
    b_with_family:      '#7E57C2',
    shared_with_family: '#78909C',
    family_only:        '#BDBDBD',
  };
  const catBorders = {
    specific_a:         '#C62828',
    specific_b:         '#6A1B9A',
    pair_exclusive:     '#00796B',
    a_with_family:      '#BF360C',
    b_with_family:      '#4527A0',
    shared_with_family: '#546E7A',
    family_only:        '#9E9E9E',
  };
  const catLabels = {
    pair_exclusive:     'Both paralogs only',
    shared_with_family: 'Both + other family',
    specific_a:         `${geneA} only`,
    specific_b:         `${geneB} only`,
    a_with_family:      `${geneA} + family (not ${geneB})`,
    b_with_family:      `${geneB} + family (not ${geneA})`,
    family_only:        'Other family only',
  };

  // Build per-row merged group indices + coverage
  const rowGroupIndices = [];
  const rowGroupCoverage = [];
  for (const row of displayRows) {
    const indices = [];
    const cov = new Map();
    for (let gi = 0; gi < nMerged; gi++) {
      const firstBi = mergedGroups[gi][0];
      let memCount = 0;
      for (const mem of row.members) { if (blocks[firstBi].positions[mem.num]) memCount++; }
      if (memCount > 0) {
        indices.push(gi);
        cov.set(gi, memCount / row.members.length);
      }
    }
    rowGroupIndices.push(indices);
    rowGroupCoverage.push(cov);
  }

  // Draw tracks
  let yPos = padTop;
  const nHumanAll = sequences.filter(s => s.is_human !== false).length;

  for (let ri = 0; ri < displayRows.length; ri++) {
    const row = displayRows[ri];
    const th = row.isPair ? pairTrackHeight : trackHeight;
    const myGroupIndices = rowGroupIndices[ri];
    const myCoverage = rowGroupCoverage[ri];

    // Label
    ctx.font = row.font;
    ctx.fillStyle = row.labelColor;
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    ctx.fillText(row.label, labelWidth - 8, yPos + th / 2);

    // Gap connectors (between merged groups only)
    const cy = yPos + th / 2;
    ctx.setLineDash([3, 3]);
    ctx.strokeStyle = row.isPair ? '#c0b69e' : '#d5d5d5';
    ctx.lineWidth = row.isPair ? 1.2 : 0.8;
    for (let k = 0; k < myGroupIndices.length - 1; k++) {
      const gi1 = myGroupIndices[k];
      const gi2 = myGroupIndices[k + 1];
      const x1 = groupColX[gi1] + groupColW[gi1];
      const x2 = groupColX[gi2];
      if (x2 > x1 + 1) { ctx.beginPath(); ctx.moveTo(x1, cy); ctx.lineTo(x2, cy); ctx.stroke(); }
    }
    ctx.setLineDash([]);

    // Draw merged groups
    for (const gi of myGroupIndices) {
      const group = mergedGroups[gi];
      const cat = blocks[group[0]].category;
      const bx = groupColX[gi];
      const bw = groupColW[gi];
      const coverage = myCoverage.get(gi) || 1;
      const firstBi = group[0];

      if (row.isGroup && coverage < 1 && row.members.length <= 4) {
        // === Branching mode: split into lanes ===
        const nMem = row.members.length;
        const laneH = Math.max(2, (th - 2) / nMem);
        for (let li = 0; li < nMem; li++) {
          const mem = row.members[li];
          const ly = yPos + 1 + li * laneH;
          if (blocks[firstBi].positions[mem.num]) {
            ctx.globalAlpha = 1.0;
            ctx.fillStyle = catColors[cat] || '#ccc';
            plmaRoundRect(ctx, bx, ly, bw, laneH - 0.5, 1);
            ctx.fill();
            ctx.strokeStyle = catBorders[cat] || '#999';
            ctx.lineWidth = 0.6;
            ctx.stroke();
          } else {
            ctx.globalAlpha = 0.15;
            ctx.fillStyle = '#888';
            plmaRoundRect(ctx, bx, ly, bw, laneH - 0.5, 1);
            ctx.fill();
            ctx.globalAlpha = 1.0;
          }
        }
      } else if (row.isGroup && coverage < 1) {
        // === Opacity mode for large groups ===
        ctx.globalAlpha = 0.3 + 0.7 * coverage;
        ctx.fillStyle = catColors[cat] || '#ccc';
        plmaRoundRect(ctx, bx, yPos + 1, bw, th - 2, 2);
        ctx.fill();
        ctx.globalAlpha = 1.0;
        ctx.setLineDash([2, 2]);
        ctx.strokeStyle = catBorders[cat] || '#999';
        ctx.lineWidth = 0.8;
        ctx.stroke();
        ctx.setLineDash([]);
      } else {
        // === Normal: singleton or full-coverage group ===
        ctx.globalAlpha = 1.0;
        ctx.fillStyle = catColors[cat] || '#ccc';
        plmaRoundRect(ctx, bx, yPos + 1, bw, th - 2, 2);
        ctx.fill();
        ctx.strokeStyle = row.isPair ? '#222' : (catBorders[cat] || '#999');
        ctx.lineWidth = row.isPair ? 1.6 : 0.8;
        ctx.stroke();
      }

      // Tooltip
      const isMerged = group.length > 1;
      let tooltipHtml;
      if (row.isGroup) {
        const memberNames = row.members.map(mkLabel);
        const present = row.members.filter(m => blocks[firstBi].positions[m.num]);
        const absent = row.members.filter(m => !blocks[firstBi].positions[m.num]);
        const pct = Math.round(coverage * 100);
        const blockLabel = isMerged
          ? `<strong>${group.length} merged blocks</strong> (${group.map(bi => blocks[bi].id).join(', ')})`
          : `<strong>${blocks[firstBi].id}</strong>`;
        tooltipHtml = `${blockLabel} · ${catLabels[cat] || cat}<br>`
          + `<span style="color:#2563eb">Group: ${memberNames.join(', ')}</span><br>`
          + `${pct}% coverage (${present.length}/${row.members.length})<br>`;
        if (absent.length > 0 && absent.length <= 5) {
          tooltipHtml += `<span style="color:#dc2626">Missing: ${absent.map(mkLabel).join(', ')}</span><br>`;
        }
        if (present.length > 0) {
          const totalAA = group.reduce((s, bi) => s + (blocks[bi].positions[present[0].num]?.length || 0), 0);
          tooltipHtml += `<span style="color:#888">${totalAA} aa total</span>`;
        }
      } else {
        const seq = row.members[0];
        if (isMerged) {
          const totalAA = group.reduce((s, bi) => s + (blocks[bi].positions[seq.num]?.length || 0), 0);
          const firstPos = blocks[group[0]].positions[seq.num];
          const lastPos = blocks[group[group.length - 1]].positions[seq.num];
          tooltipHtml = `<strong>${group.length} merged blocks</strong> (${group.map(bi => blocks[bi].id).join(', ')}) · ${catLabels[cat] || cat}<br>`
            + `${row.label}: pos ${firstPos?.start || '?'}–${lastPos?.end || '?'} (${totalAA} aa)<br>`
            + `<span style="color:#888">${blocks[firstBi].n_human || blocks[firstBi].n_seqs} of ${nHumanAll} human paralogs</span>`;
        } else {
          const pos = blocks[firstBi].positions[seq.num];
          const aaSeq = pos.seq || '';
          let aaHtml = '';
          if (aaSeq.length > 0) {
            const wrapped = aaSeq.match(/.{1,40}/g) || [aaSeq];
            aaHtml = `<code style="font-size:10px;color:#555;word-break:break-all;line-height:1.3;display:block;margin-top:3px;">${wrapped.join('<br>')}</code>`;
          }
          const nHuman = blocks[firstBi].n_human || blocks[firstBi].n_seqs;
          const nTotal = blocks[firstBi].n_seqs;
          const memberInfo = nTotal > nHuman
            ? `${nHuman} of ${nHumanAll} human paralogs (+${nTotal - nHuman} orthologs)`
            : `${nHuman} of ${nHumanAll} human paralogs`;
          tooltipHtml = `<strong>${blocks[firstBi].id}</strong> · ${catLabels[cat] || cat}<br>`
            + `${row.label}: pos ${pos.start}–${pos.end} (${pos.length} aa)<br>`
            + `<span style="color:#888">${memberInfo}</span>` + aaHtml;
        }
      }
      plmaHitRegions.push({ x: bx, y: yPos, w: bw, h: th, tooltip: tooltipHtml });
    }

    // Separator after pair tracks
    if (ri === 1 && nRows > 2) {
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
  ctx.globalAlpha = 1.0;

  // Legend
  if (legendEl) {
    const usedCats = new Set(blocks.map(b => b.category));
    let html = '';
    for (const [cat, lbl] of Object.entries(catLabels)) {
      if (!usedCats.has(cat)) continue;
      html += `<span style="display:inline-flex;align-items:center;gap:4px;">`
        + `<span style="display:inline-block;width:14px;height:10px;border-radius:2px;background:${catColors[cat]};border:1px solid ${catBorders[cat]}"></span>`
        + `<span>${lbl}</span></span>`;
    }
    html += `<span style="display:inline-flex;align-items:center;gap:4px;">`
      + `<span style="display:inline-block;width:14px;border-top:1.5px dashed #b0a890;height:0;"></span>`
      + `<span>Gap between blocks</span></span>`;
    if (!showAllParalogs && displayRows.some(r => r.isGroup)) {
      html += `<span style="display:inline-flex;align-items:center;gap:4px;">`
        + `<span style="display:inline-block;width:14px;height:10px;border-radius:2px;background:#78909C;opacity:0.45;border:1px dashed #546E7A"></span>`
        + `<span>Partial in group</span></span>`;
    }
    legendEl.innerHTML = html;
  }

  // Summary
  if (summaryEl) {
    const s = plma.summary || {};
    const nHuman = sequences.filter(x => x.is_human !== false).length;
    const nOrthoShown = showOrthologs ? orthologs.length : 0;
    const parts = [];
    if (s.shared_with_family) parts.push(`${s.shared_with_family} aa both+family`);
    if (s.pair_exclusive) parts.push(`${s.pair_exclusive} aa both only`);
    if (s.specific_a) parts.push(`${s.specific_a} aa ${geneA} only`);
    if (s.a_with_family) parts.push(`${s.a_with_family} aa ${geneA}+fam`);
    if (s.specific_b) parts.push(`${s.specific_b} aa ${geneB} only`);
    if (s.b_with_family) parts.push(`${s.b_with_family} aa ${geneB}+fam`);
    if (s.family_only) parts.push(`${s.family_only} aa family only`);

    const nMergedCount = isSimplified ? mergedGroups.filter(g => g.length > 1).length : 0;
    let seqDesc;
    if (!showAllParalogs && humanParalogs.length > 2) {
      const nGroups = displayRows.filter(r => r.isGroup).length;
      const shownHuman = displayRows.filter(r => !r.isPair).reduce((sum, r) => sum + r.members.filter(m => m.is_human !== false).length, 0);
      const nHidden = humanParalogs.length - shownHuman;
      seqDesc = `${nRows} rows (${nHuman} paralogs`;
      if (nGroups > 0) seqDesc += `, ${nGroups} grouped`;
      if (nHidden > 0) seqDesc += `, ${nHidden} hidden`;
      if (nOrthoShown > 0) seqDesc += ` + ${nOrthoShown} orthologs`;
      seqDesc += ')';
    } else {
      seqDesc = nOrthoShown > 0 ? `${nHuman} human + ${nOrthoShown} orthologs` : `${nHuman} family members`;
    }
    let blockDesc = `${blocks.length} conserved blocks`;
    if (nMergedCount > 0) blockDesc += ` (${nMergedCount} merged into ${nMerged} visual groups)`;
    summaryEl.textContent = `${blockDesc} across ${seqDesc}` +
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
  // Enable outline by default
  try {
    const pp = plugin.canvas3d?.props?.postprocessing;
    if (pp) {
      plugin.canvas3d.setProps({ postprocessing: { ...pp, outline: { name: 'on', params: { scale: 1, threshold: 0.33, color: { r: 0, g: 0, b: 0 }, includeTransparent: true } } } });
    }
  } catch(e) { console.warn('Could not enable outline:', e); }

  // Reapply color theme when Molstar representations change (e.g., from built-in controls)
  let _colorReapplyTimer = null;
  try {
    plugin.state.data.events.changed.subscribe(() => {
      if (!structureReady || !currentColorMode || currentColorMode === 'uniform') return;
      clearTimeout(_colorReapplyTimer);
      _colorReapplyTimer = setTimeout(async () => {
        try {
          const theme = themeForColorMode(currentColorMode);
          await applyColorTheme(theme);
          await renderSelections();
        } catch(e) {}
      }, 150);
    });
  } catch(e) { console.warn('Could not subscribe to Molstar state changes:', e); }
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
let hlTrackA = null, hlTrackB = null;

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
    out.push({ ...base, end:e, to:e, color, opacity, id, label, druggability: r.druggability, ss_type: r.ss_type });
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

  document.querySelectorAll('#regionsTableA tr, #regionsTableB tr, #drugTableA tr, #drugTableB tr').forEach(tr => {
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
let pendingSelectionLoci = null; // for zoom-to-selection

function focusMainViewerSelection() {
  if (!plugin || !structureReady || !pendingSelectionLoci) return;
  try { plugin.managers.camera.focusLoci(pendingSelectionLoci); } catch(e) {}
}

function setViewerLocked(locked) {
  viewerLocked = locked;
  const btn = document.getElementById('lockViewer');
  if (btn) {
    btn.innerHTML = locked ? '&#x1F512;' : '&#x1F513;';
    btn.title = locked ? 'Unlock hover highlight' : 'Lock hover highlight';
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

// Domain color palette - distinct modern colors for unique domain types
const DOMAIN_PALETTE = [
  '#4e79a7','#f28e2b','#e15759','#76b7b2','#59a14f',
  '#edc948','#b07aa1','#ff9da7','#9c755f','#bab0ac',
  '#86bcb6','#8cd17d','#b6992d','#499894','#d37295',
];
const TED_PALETTE = ['#00897b','#26a69a','#4db6ac','#80cbc4','#009688','#00796b'];

function assignDomainColors(domainsA, domainsB) {
  const nameToColor = {};
  const tedNameToColor = {};
  let domIdx = 0, tedIdx = 0;
  const allDoms = [...(domainsA||[]), ...(domainsB||[])];
  for (const d of allDoms) {
    if (d.type === 'CAV' || d.type === 'Cavity' || d.type === 'DrugCLIP') continue;
    const name = d.label || d.name || d.type || 'unknown';
    if (d.raw_type === 'TED') {
      if (!(name in tedNameToColor)) { tedNameToColor[name] = TED_PALETTE[tedIdx++ % TED_PALETTE.length]; }
      d.color = tedNameToColor[name];
    } else {
      if (!(name in nameToColor)) { nameToColor[name] = DOMAIN_PALETTE[domIdx++ % DOMAIN_PALETTE.length]; }
      d.color = nameToColor[name];
    }
  }
}

function getVisibleChains() {
  const chains = [];
  if (chainVisible.A) chains.push('A');
  if (chainVisible.B) chains.push('B');
  return chains;
}

// Filter a base64-encoded PDB to keep only the specified chain(s)
function filterPdbChains(pdb64, keepA, keepB) {
  if (keepA && keepB) return pdb64;
  const text = atob(pdb64);
  const lines = text.split('\n');
  const out = [];
  for (const line of lines) {
    if ((line.startsWith('ATOM') || line.startsWith('HETATM') || line.startsWith('TER')) && line.length > 21) {
      const chain = line[21];
      if (chain === 'A' && !keepA) continue;
      if (chain === 'B' && !keepB) continue;
    }
    out.push(line);
  }
  return btoa(out.join('\n'));
}

async function applyChainVisibility() {
  if (!plugin) return;

  let pdb64;

  // Get the base PDB (colored or full)
  if (currentColorMode && currentColorMode !== 'uniform') {
    pdb64 = getColoredPdb(currentColorMode);
  }
  if (!pdb64) {
    pdb64 = PDB64_FULL;
  }

  // Filter out hidden chain(s) at the PDB level — reliable, no Molstar hierarchy manipulation
  pdb64 = filterPdbChains(pdb64, chainVisible.A, chainVisible.B);

  await reloadViewerWith(pdb64, true);

  if (currentColorMode && currentColorMode !== 'uniform') {
    const theme = themeForColorMode(currentColorMode);
    await applyColorTheme(theme);
  }

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
      pendingSelectionLoci = null;
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

    // Build combined loci for focus/zoom (union of A and B elements)
    const allElements = [...elementsChainA, ...elementsChainB];
    if (allElements.length > 0) {
      pendingSelectionLoci = { kind: 'element-loci', structure: structureData, elements: allElements };
    } else {
      pendingSelectionLoci = null;
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
let currentPdbeColorMode = 'grey';
let pdbeComplexVisible = true;
let _pdbeResSeqToUniprot = null;
let _pdbeResSeqCacheEntry = null;

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

    // Reapply PDBe color theme when representations change
    let _pdbeColorReapplyTimer = null;
    try {
      pdbePlugin.state.data.events.changed.subscribe(() => {
        if (!pdbeStructureReady || !currentPdbeColorMode) return;
        clearTimeout(_pdbeColorReapplyTimer);
        _pdbeColorReapplyTimer = setTimeout(async () => {
          try {
            await applyPdbeCustomTheme(getPdbeColorTheme(currentPdbeColorMode));
            // Reapply both highlight (Gene B) and selection (Gene A) after state change
            if (pendingPdbeHighlightLoci || pendingPdbeSelectLoci) reapplyPdbeHighlight();
          } catch(e) {}
        }, 150);
      });
    } catch(e) { console.warn('Could not subscribe to PDBe state changes:', e); }
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

// === PDBe Coloring System ===

function extractPlddtFromPdb(b64) {
  if (!b64) return [];
  try {
    const pdb = atob(b64);
    const bfactors = {};
    for (const line of pdb.split('\n')) {
      if (line.startsWith('ATOM') && line.substring(12, 16).trim() === 'CA') {
        const resSeq = parseInt(line.substring(22, 26).trim(), 10);
        const bf = parseFloat(line.substring(60, 66).trim());
        if (!isNaN(bf) && !isNaN(resSeq)) bfactors[resSeq] = bf;
      }
    }
    const maxRes = Math.max(0, ...Object.keys(bfactors).map(Number));
    if (!maxRes) return [];
    const arr = [];
    for (let i = 1; i <= maxRes; i++) arr.push(bfactors[i] || 0);
    return arr;
  } catch(e) { return []; }
}

const AA3TO1 = {
  'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q',
  'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
  'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
  'SEC':'U','PYL':'O','MSE':'M'
};

function parseCifAtomSiteCols(cifText) {
  const atIdx = cifText.indexOf('_atom_site.');
  if (atIdx < 0) return null;
  const loopIdx = cifText.lastIndexOf('loop_', atIdx);
  const cols = {};
  let colIdx = 0;
  let pos = cifText.indexOf('\n', loopIdx) + 1;
  while (pos < cifText.length) {
    const nl = cifText.indexOf('\n', pos);
    if (nl < 0) break;
    const line = cifText.substring(pos, nl).trim();
    if (line.startsWith('_atom_site.')) {
      cols[line.substring(11)] = colIdx++;
    } else if (colIdx > 0) {
      return { cols, dataPos: pos, numCols: colIdx };
    }
    pos = nl + 1;
  }
  return null;
}

function buildPdbeResSeqMapCif(cifText, targetChainIds) {
  const info = parseCifAtomSiteCols(cifText);
  if (!info) return [];
  const { cols, dataPos } = info;
  const iGroup = cols['group_PDB'], iChain = cols['auth_asym_id'];
  const iResSeq = cols['auth_seq_id'], iCompId = cols['auth_comp_id'];
  if (iGroup === undefined || iChain === undefined || iResSeq === undefined || iCompId === undefined) return [];
  const seen = new Set(), result = [];
  let pos = dataPos;
  while (pos < cifText.length) {
    const nl = cifText.indexOf('\n', pos);
    if (nl < 0) break;
    const line = cifText.substring(pos, nl).trim();
    pos = nl + 1;
    if (!line || line.startsWith('#') || line.startsWith('_') || line.startsWith('loop_')) break;
    const fields = line.split(/\s+/);
    if (fields[iGroup] !== 'ATOM') continue;
    const chain = fields[iChain];
    if (!targetChainIds.has(chain)) continue;
    const resSeq = parseInt(fields[iResSeq], 10);
    const key = chain + '_' + resSeq;
    if (seen.has(key)) continue;
    seen.add(key);
    const aa1 = AA3TO1[fields[iCompId]];
    if (!aa1) continue;
    result.push({ resSeq, aa1, chain });
  }
  return result;
}

function buildPdbeResSeqMap(text, targetChainIds) {
  if (text.trimStart().startsWith('data_')) return buildPdbeResSeqMapCif(text, targetChainIds);
  const lines = text.split('\n'), seen = new Set(), result = [];
  for (const line of lines) {
    if (line.length < 54 || !line.startsWith('ATOM')) continue;
    const chain = line[21];
    if (!targetChainIds.has(chain)) continue;
    const resSeq = parseInt(line.substring(22, 26).trim(), 10);
    const key = chain + '_' + resSeq;
    if (seen.has(key)) continue;
    seen.add(key);
    const aa1 = AA3TO1[line.substring(17, 20).trim()];
    if (!aa1) continue;
    result.push({ resSeq, aa1, chain });
  }
  return result;
}

function mapPdbeToUniprot(resSeqArr, uniprotSeq) {
  if (!resSeqArr.length || !uniprotSeq) return {};
  const offsets = new Map();
  for (const r of resSeqArr) {
    for (let uPos = 1; uPos <= uniprotSeq.length; uPos++) {
      if (uniprotSeq[uPos - 1] === r.aa1) {
        const off = uPos - r.resSeq;
        offsets.set(off, (offsets.get(off) || 0) + 1);
      }
    }
  }
  let bestOffset = 0, bestCount = 0;
  for (const [off, cnt] of offsets) { if (cnt > bestCount) { bestCount = cnt; bestOffset = off; } }
  const mapping = {};
  for (const r of resSeqArr) {
    const uPos = r.resSeq + bestOffset;
    if (uPos >= 1 && uPos <= uniprotSeq.length) mapping[r.resSeq] = uPos;
  }
  console.log('PDBe->UniProt: offset=' + bestOffset + ', matched=' + bestCount + '/' + resSeqArr.length);
  return mapping;
}

function modifyPdbeBfactorsCif(cifText, targetChainIds, bfMap, targetDefaultBf, nonTargetBf, chainShading) {
  const info = parseCifAtomSiteCols(cifText);
  if (!info) return cifText;
  const { cols, dataPos, numCols } = info;
  const iBf = cols['B_iso_or_equiv'], iGroup = cols['group_PDB'];
  const iChain = cols['auth_asym_id'], iResSeq = cols['auth_seq_id'];
  const iCompId = cols['auth_comp_id'] !== undefined ? cols['auth_comp_id'] : cols['label_comp_id'];
  if (iBf === undefined || iGroup === undefined || iChain === undefined || iResSeq === undefined) return cifText;
  const NUCLEOTIDES = new Set(['DA','DC','DG','DT','DU','A','C','G','U','PSU','5MC','7MG','OMC','OMG','H2U','5MU','4SU','1MA','M2G']);
  let chainBf = null;
  if (chainShading && iCompId !== undefined) {
    const chainCounts = {};
    let pos = dataPos;
    while (pos < cifText.length) {
      const nl = cifText.indexOf('\n', pos); if (nl < 0) break;
      const line = cifText.substring(pos, nl).trim(); pos = nl + 1;
      if (!line || line.startsWith('#') || line.startsWith('_') || line.startsWith('loop_')) break;
      const fields = line.split(/\s+/);
      if (fields.length >= numCols && (fields[iGroup] === 'ATOM' || fields[iGroup] === 'HETATM')) {
        const chain = fields[iChain];
        if (targetChainIds.has(chain)) continue;
        if (!chainCounts[chain]) chainCounts[chain] = { nuc: 0, prot: 0 };
        if (NUCLEOTIDES.has(fields[iCompId])) chainCounts[chain].nuc++; else chainCounts[chain].prot++;
      }
    }
    const protChains = [], nucChains = [];
    for (const [ch, c] of Object.entries(chainCounts)) {
      const t = c.nuc + c.prot;
      if (t > 0 && c.nuc / t > 0.5) nucChains.push(ch); else protChains.push(ch);
    }
    chainBf = {};
    protChains.forEach((ch, i) => { chainBf[ch] = [-25, -50, -75][i % 3]; });
    nucChains.forEach(ch => { chainBf[ch] = -100; });
  }
  const header = cifText.substring(0, dataPos);
  const lines = cifText.substring(dataPos).split('\n');
  const out = [];
  for (const line of lines) {
    const trimmed = line.trim();
    if (trimmed && !trimmed.startsWith('#') && !trimmed.startsWith('_') && !trimmed.startsWith('loop_')) {
      const fields = trimmed.split(/\s+/);
      if (fields.length >= numCols && (fields[iGroup] === 'ATOM' || fields[iGroup] === 'HETATM')) {
        const chain = fields[iChain];
        if (targetChainIds.has(chain)) {
          const resSeq = parseInt(fields[iResSeq], 10);
          fields[iBf] = ((resSeq in bfMap) ? bfMap[resSeq] : targetDefaultBf).toFixed(2);
        } else if (chainBf && chainBf[chain] !== undefined) {
          fields[iBf] = chainBf[chain].toFixed(2);
        } else { fields[iBf] = nonTargetBf.toFixed(2); }
        out.push(fields.join(' ')); continue;
      }
    }
    out.push(line);
  }
  return header + out.join('\n');
}

function modifyPdbeBfactors(text, targetChainIds, bfMap, targetDefaultBf, nonTargetBf, chainShading) {
  if (targetDefaultBf === undefined) targetDefaultBf = 50;
  if (nonTargetBf === undefined) nonTargetBf = -2;
  if (chainShading === undefined) chainShading = false;
  if (text.trimStart().startsWith('data_')) return modifyPdbeBfactorsCif(text, targetChainIds, bfMap, targetDefaultBf, nonTargetBf, chainShading);
  const NUCLEOTIDES = new Set(['DA','DC','DG','DT','DU','A','C','G','U','PSU','5MC','7MG','OMC','OMG','H2U','5MU','4SU','1MA','M2G']);
  const lines = text.split('\n');
  let chainBf = null;
  if (chainShading) {
    const chainTypes = {};
    for (const line of lines) {
      if (line.length < 22 || (!line.startsWith('ATOM') && !line.startsWith('HETATM'))) continue;
      const chain = line[21];
      if (targetChainIds.has(chain)) continue;
      if (!chainTypes[chain]) chainTypes[chain] = { nuc: 0, prot: 0 };
      if (NUCLEOTIDES.has(line.substring(17, 20).trim())) chainTypes[chain].nuc++; else chainTypes[chain].prot++;
    }
    const protChains = [], nucChains = [];
    for (const [ch, c] of Object.entries(chainTypes)) {
      const t = c.nuc + c.prot;
      if (t > 0 && c.nuc / t > 0.5) nucChains.push(ch); else protChains.push(ch);
    }
    chainBf = {};
    protChains.forEach((ch, i) => { chainBf[ch] = [-25, -50, -75][i % 3]; });
    nucChains.forEach(ch => { chainBf[ch] = -100; });
  }
  const out = [];
  for (const line of lines) {
    if (line.length >= 66 && (line.startsWith('ATOM') || line.startsWith('HETATM'))) {
      const chain = line[21];
      let bf;
      if (targetChainIds.has(chain)) {
        const resSeq = parseInt(line.substring(22, 26).trim(), 10);
        bf = (resSeq in bfMap) ? bfMap[resSeq] : targetDefaultBf;
      } else if (chainBf && chainBf[chain] !== undefined) {
        bf = chainBf[chain];
      } else { bf = nonTargetBf; }
      out.push(line.substring(0, 60) + bf.toFixed(2).padStart(6) + line.substring(66));
    } else { out.push(line); }
  }
  return out.join('\n');
}

function buildPdbeGenericBfactorMap(resSeqToUniprot, mode, isGeneB) {
  const builderMap = {
    plddt: buildPlddtBfactorMaps, am: buildAmBfactorMaps, dam: buildDamBfactorMaps,
    aligned: buildAlignedBfactorMaps, domains: buildDomainBfactorMaps,
    ss: buildSsBfactorMaps, cavities: buildCavityBfactorMaps,
    drugclip: buildDrugclipBfactorMaps, plma: buildPlmaBfactorMaps,
  };
  const builder = builderMap[mode];
  if (!builder) return null;
  const maps = builder();
  if (!maps) return null;
  const geneMap = isGeneB ? maps.B : maps.A;
  if (!geneMap) return null;
  const bfMap = {};
  for (const [resSeqStr, uniprotPos] of Object.entries(resSeqToUniprot)) {
    const val = geneMap[uniprotPos];
    if (val !== undefined) {
      let sv = val;
      if (mode === 'plddt') sv = (val + 1) * 25;
      else if (mode === 'aligned' || mode === 'drugclip') sv = val * 100;
      else if (mode === 'ss') sv = val * 50;
      else if (mode === 'cavities') sv = Math.round(val * 33.33);
      else if (mode === 'plma') sv = val * 20;
      else if (mode === 'domains') { const nd = Math.max(1, Object.keys(window._domainColorNames || {}).length); sv = Math.round(val * 100 / nd); }
      bfMap[parseInt(resSeqStr)] = sv;
    }
  }
  return bfMap;
}

function getPdbeTargetDefaultBf(mode) {
  const defaults = { plddt: 25, am: 0, dam: 0, aligned: 0, domains: 0, ss: 0, cavities: 0, drugclip: 0, plma: 0 };
  return defaults[mode] !== undefined ? defaults[mode] : 0;
}

function getPdbeColorTheme(mode) {
  if (mode === 'grey') {
    // Target chain B=100 → vivid cyan; context chains B=-25,-50,-75,-100 → grey shades
    return { name: 'uncertainty', params: { domain: [-100, 100], list: { kind: 'interpolate', colors: [
      0x00ACC1, 0x26C6DA, 0x80DEEA, 0xCCEEF4,  // [0-3] B=100→25: vivid cyan → light
      0xEEEEEE,                                   // [4] B=0: neutral (rarely used)
      0xCCCCCC, 0xAAAAAA, 0x777777, 0x333333     // [5-8] B=-25→-100: light→dark grey
    ]}}};
  }
  if (mode === 'plddt') {
    return { name: 'uncertainty', params: { domain: [-100, 100], list: { kind: 'interpolate', colors: [
      0x0053d6, 0x65cbf3, 0xffdb13, 0xff7d45, 0xff7d45,
      0xD4C5A9, 0xC0A882, 0x8B7355, 0x1a1a1a
    ]}}};
  }
  if (mode === 'am' || mode === 'dam') {
    return { name: 'uncertainty', params: { domain: [-100, 100], list: { kind: 'interpolate', colors: [
      0xd62728, 0xff7d45, 0xbbbbbb, 0xdddddd, 0xdddddd,
      0xD4C5A9, 0xC0A882, 0x8B7355, 0x1a1a1a
    ]}}};
  }
  if (mode === 'aligned') {
    // 4-way: 0=gap(burgundy,B=0), 25=identical(grey,B=25), 50=conservative(teal,B=50), 100=radical(green,B=100)
    // Reversed: high B(100)→first color(index 0), low B(-100)→last color(index 8)
    return { name: 'uncertainty', params: { domain: [-100, 100], list: { kind: 'interpolate', colors: [
      0x4CAF50, 0x26A69A, 0x009688, 0xAAAAAA, 0x8B1A1A,
      0xD4C5A9, 0xC0A882, 0x8B7355, 0x1a1a1a
    ]}}};
  }
  if (mode === 'domains') {
    // data: 0=no-domain(grey), scaled to [100/N, 200/N, ..., 100]; chain shading at -25,-50,-75,-100
    return { name: 'uncertainty', params: { domain: [-100, 100], list: { kind: 'interpolate', colors: [
      0x1565C0, 0x1976D2, 0x42A5F5, 0x90CAF9, 0xcccccc,
      0xD4C5A9, 0xC0A882, 0x8B7355, 0x1a1a1a
    ]}}};
  }
  if (mode === 'ss') {
    // data: 0=coil(grey), 50=strand(yellow), 100=helix(red); chain shading at -25,-50,-75,-100
    return { name: 'uncertainty', params: { domain: [-100, 100], list: { kind: 'interpolate', colors: [
      0xFF0066, 0xE06020, 0xFFCC00, 0xCCCC88, 0xdddddd,
      0xD4C5A9, 0xC0A882, 0x8B7355, 0x1a1a1a
    ]}}};
  }
  if (mode === 'cavities') {
    // data: 0=none(grey), 33=weak, 67=medium, 100=strong; chain shading at -25,-50,-75,-100
    return { name: 'uncertainty', params: { domain: [-100, 100], list: { kind: 'interpolate', colors: [
      0xe65100, 0xff6d00, 0xff9800, 0xffcc80, 0xdddddd,
      0xD4C5A9, 0xC0A882, 0x8B7355, 0x1a1a1a
    ]}}};
  }
  if (mode === 'drugclip') {
    // data: 0=no-pocket(white), 100=pocket(red); chain shading at -25,-50,-75,-100
    return { name: 'uncertainty', params: { domain: [-100, 100], list: { kind: 'interpolate', colors: [
      0xc62828, 0xc62828, 0xc62828, 0xEE9999, 0xf7f7f7,
      0xD4C5A9, 0xC0A882, 0x8B7355, 0x1a1a1a
    ]}}};
  }
  if (mode === 'plma') {
    // data: 0=none, 20=family-only, ..., 100=specific; chain shading at -25,-50,-75,-100
    return { name: 'uncertainty', params: { domain: [-100, 100], list: { kind: 'interpolate', colors: [
      0xEF5350, 0xFFA726, 0x26A69A, 0xFFCA28, 0xBDBDBD,
      0xD4C5A9, 0xC0A882, 0x8B7355, 0x1a1a1a
    ]}}};
  }
  return { name: 'uniform', params: { value: 0xcccccc } };
}

function getPdbeColorLegendHtml(mode) {
  const sw = (c, lbl) => '<span style="display:inline-flex;align-items:center;gap:3px;margin-right:5px"><span style="display:inline-block;width:11px;height:11px;background:' + c + ';border:1px solid #ccc;border-radius:2px;flex-shrink:0"></span>' + lbl + '</span>';
  const legends = {
    plddt:    '<strong>pLDDT:</strong> ' + sw('#0053d6','&gt;90') + sw('#65cbf3','70-90') + sw('#ffdb13','50-70') + sw('#ff7d45','&le;50'),
    am:       '<strong>AlphaMissense:</strong> ' + sw('#d62728','Pathogenic') + sw('#ff7d45','Ambiguous') + sw('#bbbbbb','Benign'),
    dam:      '<strong>&Delta; AlphaMissense:</strong> ' + sw('#d62728','High') + sw('#ff7d45','Med') + sw('#bbbbbb','Low'),
    aligned:  '<strong>Substitution:</strong> ' + sw('#4CAF50','Radical') + sw('#009688','Conservative') + sw('#AAAAAA','Identical') + sw('#8B1A1A','Gap'),
    domains:  '<strong>Domains:</strong> colored by type',
    ss:       '<strong>2D Structure:</strong> ' + sw('#FF0066','&alpha;-helix') + sw('#FFCC00','&beta;-strand') + sw('#dddddd','Coil'),
    cavities: '<strong>Cavities:</strong> ' + sw('#e65100','Strong') + sw('#ff9800','Medium') + sw('#ffc107','Weak'),
    drugclip: '<strong>DrugCLIP:</strong> ' + sw('#c62828','Pocket') + sw('#f7f7f7','None'),
    plma:     '<strong>PLMA:</strong> ' + sw('#EF5350','Specific') + sw('#FFA726','+Family') + sw('#26A69A','Pair-excl') + sw('#FFCA28','Shared') + sw('#BDBDBD','Family-only'),
  };
  const chainNote = sw('#C9B99A','Context') + sw('#1a1a1a','DNA/RNA');
  return (legends[mode] || '') + ' | ' + chainNote;
}

async function applyPdbeCustomTheme(theme) {
  if (!pdbePlugin || !pdbeStructureReady) return;
  const snapshot = pdbePlugin.canvas3d ? pdbePlugin.canvas3d.camera.getSnapshot() : null;
  try {
    const structures = pdbePlugin.managers.structure.hierarchy.current.structures;
    if (!structures || !structures.length) return;
    const update = pdbePlugin.state.data.build();
    for (const struct of structures) {
      for (const comp of (struct.components || [])) {
        for (const repr of (comp.representations || [])) {
          update.to(repr.cell).update(function(old) { return Object.assign({}, old, { colorTheme: theme }); });
        }
      }
    }
    await update.commit();
  } catch(e) { console.error('applyPdbeCustomTheme:', e); }
  if (snapshot && pdbePlugin && pdbePlugin.canvas3d && pdbePlugin.canvas3d.camera) pdbePlugin.canvas3d.camera.setState(snapshot);
}

async function zoomToPdbeProtein() {
  if (!pdbePlugin || !pdbeStructureReady || !currentPdbeEntry) return;
  const targetChains = getTargetChainIds(currentPdbeEntry);
  if (!targetChains.size) return;
  try {
    const structures = pdbePlugin.managers.structure.hierarchy.current.structures;
    if (!structures || !structures.length) return;
    const structureData = structures[0].cell?.obj?.data;
    if (!structureData) { try { pdbeViewer.resetCamera(); } catch(e) {} return; }

    const units = structureData.units || [];
    const matchingElements = [];
    for (let ui = 0; ui < units.length; ui++) {
      const unit = units[ui];
      const model = unit.model;
      if (!model) continue;
      const chains = model.atomicHierarchy?.chains;
      const chainAtomSegments = model.atomicHierarchy?.chainAtomSegments;
      if (!chains || !chainAtomSegments) continue;
      const unitElements = unit.elements;
      if (!unitElements || !unitElements.length) continue;
      const firstAtomIdx = unitElements[0];
      const chainOffsets = chainAtomSegments.offsets;
      let unitChainId = null;
      for (let ci = 0; ci < chainOffsets.length - 1; ci++) {
        if (firstAtomIdx >= chainOffsets[ci] && firstAtomIdx < chainOffsets[ci + 1]) {
          unitChainId = chains.auth_asym_id?.value(ci);
          break;
        }
      }
      if (unitChainId && targetChains.has(unitChainId)) {
        const indices = [];
        for (let i = 0; i < unitElements.length; i++) indices.push(i);
        matchingElements.push({ unit, indices });
      }
    }
    if (matchingElements.length > 0) {
      const loci = { kind: 'element-loci', structure: structureData, elements: matchingElements };
      pdbePlugin.managers.camera.focusLoci(loci);
    } else {
      try { pdbeViewer.resetCamera(); } catch(e) {}
    }
  } catch(e) {
    console.warn('zoomToPdbeProtein error:', e);
    try { pdbeViewer.resetCamera(); } catch(e2) {}
  }
}

// Pending PDBe highlight loci for hover-lock reapplication and focus
let pendingPdbeHighlightLoci = null;  // Gene B (pink highlight)
let pendingPdbeSelectLoci = null;     // Gene A (green select)
let pendingPdbeFocusLoci = null;      // for pdbeFocusSelection button
let pdbeHoverLockInterval = null;

function focusPdbeSelection() {
  if (!pdbePlugin || !pdbeStructureReady || !pendingPdbeFocusLoci) return;
  try { pdbePlugin.managers.camera.focusLoci(pendingPdbeFocusLoci); } catch(e) {}
}

function reapplyPdbeHighlight() {
  if (!pdbePlugin) return;
  if (pendingPdbeHighlightLoci) {
    try { pdbePlugin.managers.interactivity.lociHighlights.highlight({ loci: pendingPdbeHighlightLoci }); } catch(e) {}
  }
  if (pendingPdbeSelectLoci) {
    try { pdbePlugin.managers.interactivity.lociSelects.select({ loci: pendingPdbeSelectLoci }); } catch(e) {}
  }
}

function setupPdbeHoverLock() {
  if (pdbeHoverLockInterval) return;
  pdbeHoverLockInterval = setInterval(() => { reapplyPdbeHighlight(); }, 80);
}

function syncSelectionsToPdbe() {
  if (!pdbePlugin || !pdbeStructureReady || !currentPdbeEntry || !_pdbeResSeqToUniprot) return;
  const sourceAcc = currentPdbeEntry.source_acc || currentPdbeEntry.sourceAcc || '';
  const isGeneB = !!(SUMMARY && SUMMARY.gene2 && SUMMARY.gene2.uniprot === sourceAcc);
  const targetChain = isGeneB ? chainIdB : chainIdA;
  const selPositions = new Set();
  for (const entry of selection.values()) {
    if (entry.chain === targetChain) {
      for (let p = entry.start; p <= entry.end; p++) selPositions.add(p);
    }
  }
  try {
    // Clear both select and highlight first so selection always replaces default
    pdbePlugin.managers.interactivity.lociSelects.deselectAll();
    pdbePlugin.managers.interactivity.lociHighlights.clearHighlights();
    pendingPdbeHighlightLoci = null;
    pendingPdbeSelectLoci = null;
    if (!selPositions.size) return;
    const targetAuthSeqIds = new Set();
    for (const [rs, up] of Object.entries(_pdbeResSeqToUniprot)) {
      if (selPositions.has(up)) targetAuthSeqIds.add(parseInt(rs));
    }
    if (!targetAuthSeqIds.size) return;
    const structures = pdbePlugin.managers.structure.hierarchy.current.structures;
    if (!structures?.length) return;
    const sd = structures[0].cell?.obj?.data;
    if (!sd) return;
    const targetChains = getTargetChainIds(currentPdbeEntry);
    const matchingElements = [];
    for (const unit of (sd.units || [])) {
      const ah = unit.model?.atomicHierarchy;
      if (!ah) continue;
      const ue = unit.elements;
      if (!ue?.length) continue;
      const chainOffsets = ah.chainAtomSegments?.offsets;
      const chains = ah.chains;
      if (!chainOffsets || !chains) continue;
      const firstAtom = ue[0];
      let unitChainId = null;
      for (let ci = 0; ci < chainOffsets.length - 1; ci++) {
        if (firstAtom >= chainOffsets[ci] && firstAtom < chainOffsets[ci + 1]) {
          unitChainId = chains.auth_asym_id?.value(ci); break;
        }
      }
      if (!unitChainId || !targetChains.has(unitChainId)) continue;
      const resOffsets = ah.residueAtomSegments?.offsets;
      const authSeqIds = ah.residues?.auth_seq_id;
      if (!resOffsets || !authSeqIds) continue;
      const gToPos = new Map();
      for (let i = 0; i < ue.length; i++) gToPos.set(ue[i], i);
      const indices = [];
      for (let ri = 0, nRes = resOffsets.length - 1; ri < nRes; ri++) {
        if (!targetAuthSeqIds.has(authSeqIds.value(ri))) continue;
        for (let a = resOffsets[ri]; a < resOffsets[ri + 1]; a++) {
          const pos = gToPos.get(a); if (pos !== undefined) indices.push(pos);
        }
      }
      if (indices.length) matchingElements.push({ unit, indices });
    }
    if (matchingElements.length) {
      const loci = { kind: 'element-loci', structure: sd, elements: matchingElements };
      pendingPdbeFocusLoci = loci;  // save for Focus button
      if (isGeneB) {
        // Gene B: pink highlight — lock via interval to prevent fading on hover
        pdbePlugin.managers.interactivity.lociHighlights.highlight({ loci });
        pendingPdbeHighlightLoci = loci;
        pendingPdbeSelectLoci = null;
        setupPdbeHoverLock();
      } else {
        // Gene A: green select (more persistent than highlight)
        pdbePlugin.managers.interactivity.lociSelects.select({ loci });
        pendingPdbeSelectLoci = loci;
        pendingPdbeHighlightLoci = null;
      }
    } else {
      pendingPdbeFocusLoci = null;
    }
  } catch(e) { console.warn('syncSelectionsToPdbe error:', e); }
}

async function togglePdbeComplexVisibility(visible) {
  if (!pdbePlugin || !pdbeStructureReady) return;
  const snapshot = pdbePlugin.canvas3d ? pdbePlugin.canvas3d.camera.getSnapshot() : null;
  try {
    const structures = pdbePlugin.managers.structure.hierarchy.current.structures;
    if (!structures || !structures.length) return;
    const pdbeStruct = structures[0];
    for (const comp of (pdbeStruct.components || [])) {
      for (const repr of (comp.representations || [])) {
        pdbePlugin.state.data.updateCellState(repr.cell.transform.ref, { isHidden: !visible });
      }
    }
    pdbeComplexVisible = visible;
  } catch(e) { console.error('togglePdbeComplexVisibility:', e); }
  if (snapshot && pdbePlugin && pdbePlugin.canvas3d && pdbePlugin.canvas3d.camera) pdbePlugin.canvas3d.camera.setState(snapshot);
}

async function applyPdbeColorMode(mode) {
  if (!currentPdbeEntry) return;
  currentPdbeColorMode = mode;
  const entry = currentPdbeEntry;
  const b64 = entry.coord_b64 || entry.coordB64 || entry.pdb_b64 || entry.pdbB64 || '';
  const format = entry.coord_format || entry.coordFormat || 'pdb';
  const pdbId = (entry.pdb_id || entry.pdbId || '').toLowerCase();
  const geneName = getGeneNameForAccession(entry.source_acc || entry.sourceAcc || '');

  let cameraSnapshot = null;
  if (pdbePlugin && pdbePlugin.canvas3d && pdbePlugin.canvas3d.camera) {
    try { cameraSnapshot = pdbePlugin.canvas3d.camera.getSnapshot(); } catch(e) {}
  }

  if (!mode || mode === 'grey') {
    // Color protein of interest in vivid cyan, context chains in grey shades
    if (b64 && b64.length > 100) {
      const pdbText = atob(b64);
      const targetChains = getTargetChainIds(entry);
      // Target chain: B=100 (vivid), non-target: B=-25..-100 (grey shading per chain)
      const modifiedText = modifyPdbeBfactors(pdbText, targetChains, {}, 100, -25, true);
      try { await pdbePlugin?.clear(); } catch(e) {}
      pdbeStructureReady = false;
      const isCif = modifiedText.trimStart().startsWith('data_');
      const blob = new Blob([modifiedText], { type: isCif ? 'chemical/x-mmcif' : 'chemical/x-pdb' });
      const url = URL.createObjectURL(blob);
      await pdbeViewer.loadStructureFromUrl(url, isCif ? 'mmcif' : 'pdb');
      pdbeStructureReady = true;
      URL.revokeObjectURL(url);
      await new Promise(r => setTimeout(r, 200));
      await applyPdbeCustomTheme(getPdbeColorTheme('grey'));
    } else if (pdbId) {
      await fetchAndLoadPdbeStructure(pdbId);
    }
    if (cameraSnapshot && pdbePlugin?.canvas3d?.camera) {
      try { pdbePlugin.canvas3d.camera.setState(cameraSnapshot); } catch(e) {}
    } else {
      // First load — zoom to the target protein chain so it's centered and clipping is right
      await zoomToPdbeProtein();
    }
    updatePdbeLegend('Context chains in <strong>grey</strong>. <strong>' + (geneName || 'Protein') + '</strong> highlighted in <span style="color:#00ACC1">teal</span>. Click "Highlight Protein" to highlight specific region.');
    return;
  }

  if (!b64 || b64.length < 10) {
    updatePdbeLegend('<strong>Error:</strong> No structure data available for coloring.');
    return;
  }

  updatePdbeLegend('Applying coloring...');
  await initPdbeMolstar();
  if (!pdbeViewer || !pdbePlugin) return;

  const pdbText = atob(b64);
  const targetChains = getTargetChainIds(entry);

  if (!_pdbeResSeqToUniprot || _pdbeResSeqCacheEntry !== entry) {
    const sourceAcc = entry.source_acc || entry.sourceAcc || '';
    const isGeneB2 = !!(SUMMARY && SUMMARY.gene2 && SUMMARY.gene2.uniprot === sourceAcc);
    const uniprotSeq = isGeneB2 ? (DATA.taln || '').replace(/-/g, '') : (DATA.qaln || '').replace(/-/g, '');
    _pdbeResSeqToUniprot = mapPdbeToUniprot(buildPdbeResSeqMap(pdbText, targetChains), uniprotSeq);
    _pdbeResSeqCacheEntry = entry;
  }

  const sourceAcc2 = entry.source_acc || entry.sourceAcc || '';
  const isGeneB = !!(SUMMARY && SUMMARY.gene2 && SUMMARY.gene2.uniprot === sourceAcc2);
  const bfMap = buildPdbeGenericBfactorMap(_pdbeResSeqToUniprot, mode, isGeneB);

  if (!bfMap || !Object.keys(bfMap).length) {
    updatePdbeLegend('<strong>Note:</strong> No data available for this coloring mode.');
    return;
  }

  const targetDefaultBf = getPdbeTargetDefaultBf(mode);
  const modifiedText = modifyPdbeBfactors(pdbText, targetChains, bfMap, targetDefaultBf, -25, true);

  try { await pdbePlugin.clear(); } catch(e) {}
  pdbeStructureReady = false;
  const isCif = modifiedText.trimStart().startsWith('data_');
  const blob = new Blob([modifiedText], { type: isCif ? 'chemical/x-mmcif' : 'chemical/x-pdb' });
  const url = URL.createObjectURL(blob);
  await pdbeViewer.loadStructureFromUrl(url, isCif ? 'mmcif' : 'pdb');
  pdbeStructureReady = true;
  URL.revokeObjectURL(url);

  await new Promise(function(r) { setTimeout(r, 300); });

  await applyPdbeCustomTheme(getPdbeColorTheme(mode));
  updatePdbeLegend(getPdbeColorLegendHtml(mode));

  if (cameraSnapshot && pdbePlugin && pdbePlugin.canvas3d && pdbePlugin.canvas3d.camera) {
    try { pdbePlugin.canvas3d.camera.setState(cameraSnapshot); } catch(e) {}
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
    const descEl = document.getElementById('pdbeDescription');
    if (descEl) descEl.textContent = '-';
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
  const descEl = document.getElementById('pdbeDescription');
  if (descEl) descEl.textContent = entry.title || '-';
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
  // Reset color mode and residue cache when switching structures
  currentPdbeColorMode = 'grey';
  _pdbeResSeqToUniprot = null;
  _pdbeResSeqCacheEntry = null;
  pdbeComplexVisible = true;
  const colorSel = document.getElementById('pdbeColorBy');
  if (colorSel) colorSel.value = 'grey';
  const visChk = document.getElementById('togglePdbeComplex');
  if (visChk) visChk.checked = true;

  // Update button text with gene name
  const sourceAcc = entry.source_acc || entry.sourceAcc || '';
  const geneName = getGeneNameForAccession(sourceAcc);
  const proteinNameSpan = document.getElementById('pdbeProteinName');
  if (proteinNameSpan) {
    proteinNameSpan.textContent = geneName || 'Protein';
  }

  updatePdbeInfoBox(entry);
  updatePdbeLegend(`Loading structure...`);
  updateHighlightButtonState(false);

  // Use applyPdbeColorMode('grey') which now colors the target protein in vivid teal
  try {
    await initPdbeMolstar();
    await applyPdbeColorMode('grey');
    // Sync any existing main viewer selections to PDBe immediately after load
    try { syncSelectionsToPdbe(); } catch(e2) {}
  } catch(e) {
    console.error('Failed to load PDBe structure for entry:', entry, e);
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
    const geneName = getGeneNameForAccession(sourceAcc);

    let label = pdbId;
    if (geneName) label += ` – ${geneName}`;
    else if (sourceAcc) label += ` [${sourceAcc}]`;
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

  const zoomProteinBtn = document.getElementById('pdbeZoomProtein');
  if (zoomProteinBtn) {
    zoomProteinBtn.addEventListener('click', async () => {
      await zoomToPdbeProtein();
    }, { passive: true });
  }

  const pdbeFocusBtn = document.getElementById('pdbeFocusSelection');
  if (pdbeFocusBtn) {
    pdbeFocusBtn.addEventListener('click', () => { focusPdbeSelection(); }, { passive: true });
  }

  // PDBe color-by dropdown
  const pdbeColorBySelect = document.getElementById('pdbeColorBy');
  if (pdbeColorBySelect) {
    pdbeColorBySelect.addEventListener('change', async (e) => {
      await applyPdbeColorMode(e.target.value);
    }, { passive: true });
  }

  // PDBe complex visibility toggle
  const pdbeVisToggle = document.getElementById('togglePdbeComplex');
  if (pdbeVisToggle) {
    pdbeVisToggle.addEventListener('change', async (e) => {
      await togglePdbeComplexVisibility(e.target.checked);
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
  updateNightingaleHighlights();
  await applyMolstarSelection();
  syncSelectionsToPdbe();
}

function updateNightingaleHighlights() {
  if (!hlTrackA || !hlTrackB || !DATA) return;
  const sel = getAllSelections();
  const qmap = DATA.qposToCol || {};
  const tmap = DATA.tposToCol || {};
  const alnLen = Math.max(1, (DATA.qaln || '').length);

  function buildFeatures(chain, posToCol, color) {
    const ranges = sel.filter(s => s.chain === chain);
    if (!ranges.length) return [];
    const features = [];
    for (const r of ranges) {
      let minCol = Infinity, maxCol = -Infinity;
      for (let p = r.start; p <= r.end; p++) {
        const col = posToCol[p];
        if (col) { minCol = Math.min(minCol, col); maxCol = Math.max(maxCol, col); }
      }
      if (minCol <= maxCol) {
        features.push({ start: minCol, end: maxCol, color, opacity: 0.8 });
      }
    }
    return sanitizeRects(features, alnLen);
  }

  const featA = buildFeatures(chainIdA, qmap, '#00e676');
  const featB = buildFeatures(chainIdB, tmap, '#ff4081');

  // Force Nightingale re-render: set data then call requestUpdate if available
  hlTrackA.data = featA;
  hlTrackB.data = featB;
  if (typeof hlTrackA.requestUpdate === 'function') hlTrackA.requestUpdate();
  if (typeof hlTrackB.requestUpdate === 'function') hlTrackB.requestUpdate();

}

function toggleFeature(dom, chain) {
  if (!dom || !dom.uid) return;
  const key = selectionKey(chain, dom.uid);
  if (selection.has(key)) {
    selection.delete(key);
  } else {
    const color = dom.color
      || ((dom.type === 'DrugCLIP') ? '#c62828'
        : (dom.type === 'Cavity') ? '#ff9800'
        : (dom.raw_type === 'TED') ? '#00897b'
        : '#ffdb13');
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
  syncHighlightDropdowns();
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

  // Build position-to-column maps if not yet set
  if (!DATA.qposToCol || !Object.keys(DATA.qposToCol).length) {
    const qposToCol = {}, tposToCol = {};
    let qp = 0, tp = 0;
    for (let col = 0; col < alnLen; col++) {
      if ((DATA.qaln||'')[col] !== '-') { qp++; qposToCol[qp] = col + 1; }
      if ((DATA.taln||'')[col] !== '-') { tp++; tposToCol[tp] = col + 1; }
    }
    DATA.qposToCol = qposToCol;
    DATA.tposToCol = tposToCol;
  }

  amMatrixTracksA = [];
  amMatrixTracksB = [];

  const nav = document.createElement('nightingale-navigation');
  addRow(tbl, '', nav, 40);

  // === Sel (Highlights) ===
  hlTrackA = document.createElement('nightingale-track');
  addRow(tbl, 'Sel. '+DATA.g1, hlTrackA, 8);
  hlTrackB = document.createElement('nightingale-track');
  addRow(tbl, 'Sel. '+DATA.g2, hlTrackB, 8);

  // === Sequences ===
  const seqA = document.createElement('nightingale-sequence');
  addRow(tbl, DATA.g1+' ('+DATA.a1+') aligned', seqA, 28);
  const seqB = document.createElement('nightingale-sequence');
  addRow(tbl, DATA.g2+' ('+DATA.a2+') aligned', seqB, 28);

  // === Domains ===
  const domA = document.createElement('nightingale-track');
  addRow(tbl, 'Domains '+DATA.g1, domA, 16); trackRefs['domA'] = domA;
  const domB = document.createElement('nightingale-track');
  addRow(tbl, 'Domains '+DATA.g2, domB, 16); trackRefs['domB'] = domB;

  // === TED ===
  const tedA = document.createElement('nightingale-track');
  const tedARow = addRow(tbl, 'TED '+DATA.g1, tedA, 16); trackRefs['tedA'] = tedA;
  const tedB = document.createElement('nightingale-track');
  const tedBRow = addRow(tbl, 'TED '+DATA.g2, tedB, 16); trackRefs['tedB'] = tedB;

  // === PLMA A (expandable overview + sub-tracks) ===
  const plmaCatA = document.createElement('nightingale-track');
  const plmaARow = addRow(tbl, 'PLMA '+DATA.g1, plmaCatA, 14); trackRefs['plmaCatA'] = plmaCatA;
  const plmaSubRowsA = [];
  const plmaToggleA = document.createElement('button');
  plmaToggleA.textContent = '▼'; plmaToggleA.className = 'am-matrix-toggle';
  plmaARow.labelCell.appendChild(plmaToggleA);

  const plmaSharedA = document.createElement('nightingale-track');
  const psARow = addRow(tbl, 'Both + family', plmaSharedA, 10);
  psARow.row.style.display = 'none'; plmaSubRowsA.push(psARow.row); trackRefs['plmaSharedA'] = plmaSharedA;
  const plmaPairA = document.createElement('nightingale-track');
  const ppARow = addRow(tbl, 'Both only', plmaPairA, 10);
  ppARow.row.style.display = 'none'; plmaSubRowsA.push(ppARow.row); trackRefs['plmaPairA'] = plmaPairA;
  const plmaSpecA = document.createElement('nightingale-track');
  const pspARow = addRow(tbl, DATA.g1+' only', plmaSpecA, 10);
  pspARow.row.style.display = 'none'; plmaSubRowsA.push(pspARow.row); trackRefs['plmaSpecA'] = plmaSpecA;
  const plmaFamA = document.createElement('nightingale-track');
  const pfARow = addRow(tbl, DATA.g1+' + family', plmaFamA, 10);
  pfARow.row.style.display = 'none'; plmaSubRowsA.push(pfARow.row); trackRefs['plmaFamA'] = plmaFamA;

  plmaToggleA.addEventListener('click', () => {
    const hidden = plmaSubRowsA[0].style.display === 'none';
    plmaSubRowsA.forEach(r => { r.style.display = hidden ? '' : 'none'; });
    plmaToggleA.textContent = hidden ? '▲' : '▼';
  }, { passive: true });

  // === PLMA B (expandable overview + sub-tracks) ===
  const plmaCatB = document.createElement('nightingale-track');
  const plmaBRow = addRow(tbl, 'PLMA '+DATA.g2, plmaCatB, 14); trackRefs['plmaCatB'] = plmaCatB;
  const plmaSubRowsB = [];
  const plmaToggleB = document.createElement('button');
  plmaToggleB.textContent = '▼'; plmaToggleB.className = 'am-matrix-toggle';
  plmaBRow.labelCell.appendChild(plmaToggleB);

  const plmaSharedB = document.createElement('nightingale-track');
  const psBRow = addRow(tbl, 'Both + family', plmaSharedB, 10);
  psBRow.row.style.display = 'none'; plmaSubRowsB.push(psBRow.row); trackRefs['plmaSharedB'] = plmaSharedB;
  const plmaPairB = document.createElement('nightingale-track');
  const ppBRow = addRow(tbl, 'Both only', plmaPairB, 10);
  ppBRow.row.style.display = 'none'; plmaSubRowsB.push(ppBRow.row); trackRefs['plmaPairB'] = plmaPairB;
  const plmaSpecB = document.createElement('nightingale-track');
  const pspBRow = addRow(tbl, DATA.g2+' only', plmaSpecB, 10);
  pspBRow.row.style.display = 'none'; plmaSubRowsB.push(pspBRow.row); trackRefs['plmaSpecB'] = plmaSpecB;
  const plmaFamB = document.createElement('nightingale-track');
  const pfBRow = addRow(tbl, DATA.g2+' + family', plmaFamB, 10);
  pfBRow.row.style.display = 'none'; plmaSubRowsB.push(pfBRow.row); trackRefs['plmaFamB'] = plmaFamB;

  plmaToggleB.addEventListener('click', () => {
    const hidden = plmaSubRowsB[0].style.display === 'none';
    plmaSubRowsB.forEach(r => { r.style.display = hidden ? '' : 'none'; });
    plmaToggleB.textContent = hidden ? '▲' : '▼';
  }, { passive: true });

  // === AlphaMissense A ===
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

  // === ΔAM ===
  const dam = document.createElement('nightingale-track');
  const damRow = addRow(tbl, 'ΔAM', dam, 18);
  damRow.row.classList.add('am-main-row');
  damRow.labelCell.classList.add('am-main-label');
  damTrack = dam;

  // === AlphaMissense B ===
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

  // === Disordered ===
  const disorderA = document.createElement('nightingale-track');
  const disorderARow = addRow(tbl, 'Disordered '+DATA.g1, disorderA, 16); trackRefs['disorderA'] = disorderA;
  const disorderB = document.createElement('nightingale-track');
  const disorderBRow = addRow(tbl, 'Disordered '+DATA.g2, disorderB, 16); trackRefs['disorderB'] = disorderB;

  // === Alignment substitution type (gap/conservative/radical) ===
  const alndA = document.createElement('nightingale-track');
  const alndARow = addRow(tbl, 'Aln '+DATA.g1, alndA, 10); trackRefs['alndA'] = alndA;
  const alndB = document.createElement('nightingale-track');
  const alndBRow = addRow(tbl, 'Aln '+DATA.g2, alndB, 10); trackRefs['alndB'] = alndB;

  // === Alpha helix ===
  const helixA = document.createElement('nightingale-track');
  const helixARow = addRow(tbl, 'α '+DATA.g1, helixA, 12); trackRefs['helixA'] = helixA;
  const helixB = document.createElement('nightingale-track');
  const helixBRow = addRow(tbl, 'α '+DATA.g2, helixB, 12); trackRefs['helixB'] = helixB;

  // === Beta strand ===
  const strandA = document.createElement('nightingale-track');
  const strandARow = addRow(tbl, 'β '+DATA.g1, strandA, 12); trackRefs['strandA'] = strandA;
  const strandB = document.createElement('nightingale-track');
  const strandBRow = addRow(tbl, 'β '+DATA.g2, strandB, 12); trackRefs['strandB'] = strandB;

  // === Cavities A (expandable) ===
  const cavA = document.createElement('nightingale-track');
  const cavARow = addRow(tbl, 'Cavities '+DATA.g1, cavA, 16); trackRefs['cavA'] = cavA;
  const cavSubRowsA = [];
  const cavToggleA = document.createElement('button');
  cavToggleA.textContent = '▼'; cavToggleA.className = 'am-matrix-toggle';
  cavARow.labelCell.appendChild(cavToggleA);

  const cavStrongA = document.createElement('nightingale-track');
  const csARow = addRow(tbl, 'Strong', cavStrongA, 10);
  csARow.row.style.display = 'none'; cavSubRowsA.push(csARow.row); trackRefs['cavStrongA'] = cavStrongA;
  const cavMediumA = document.createElement('nightingale-track');
  const cmARow = addRow(tbl, 'Medium', cavMediumA, 10);
  cmARow.row.style.display = 'none'; cavSubRowsA.push(cmARow.row); trackRefs['cavMediumA'] = cavMediumA;
  const cavWeakA = document.createElement('nightingale-track');
  const cwARow = addRow(tbl, 'Weak', cavWeakA, 10);
  cwARow.row.style.display = 'none'; cavSubRowsA.push(cwARow.row); trackRefs['cavWeakA'] = cavWeakA;

  cavToggleA.addEventListener('click', () => {
    const hidden = cavSubRowsA[0].style.display === 'none';
    cavSubRowsA.forEach(r => { r.style.display = hidden ? '' : 'none'; });
    cavToggleA.textContent = hidden ? '▲' : '▼';
  }, { passive: true });

  // === Cavities B (expandable) ===
  const cavB = document.createElement('nightingale-track');
  const cavBRow = addRow(tbl, 'Cavities '+DATA.g2, cavB, 16); trackRefs['cavB'] = cavB;
  const cavSubRowsB = [];
  const cavToggleB = document.createElement('button');
  cavToggleB.textContent = '▼'; cavToggleB.className = 'am-matrix-toggle';
  cavBRow.labelCell.appendChild(cavToggleB);

  const cavStrongB = document.createElement('nightingale-track');
  const csBRow = addRow(tbl, 'Strong', cavStrongB, 10);
  csBRow.row.style.display = 'none'; cavSubRowsB.push(csBRow.row); trackRefs['cavStrongB'] = cavStrongB;
  const cavMediumB = document.createElement('nightingale-track');
  const cmBRow = addRow(tbl, 'Medium', cavMediumB, 10);
  cmBRow.row.style.display = 'none'; cavSubRowsB.push(cmBRow.row); trackRefs['cavMediumB'] = cavMediumB;
  const cavWeakB = document.createElement('nightingale-track');
  const cwBRow = addRow(tbl, 'Weak', cavWeakB, 10);
  cwBRow.row.style.display = 'none'; cavSubRowsB.push(cwBRow.row); trackRefs['cavWeakB'] = cavWeakB;

  cavToggleB.addEventListener('click', () => {
    const hidden = cavSubRowsB[0].style.display === 'none';
    cavSubRowsB.forEach(r => { r.style.display = hidden ? '' : 'none'; });
    cavToggleB.textContent = hidden ? '▲' : '▼';
  }, { passive: true });

  // === DrugCLIP ===
  const dcA = document.createElement('nightingale-track');
  const dcARow = addRow(tbl, 'DrugCLIP '+DATA.g1, dcA, 16); trackRefs['dcA'] = dcA;
  const dcB = document.createElement('nightingale-track');
  const dcBRow = addRow(tbl, 'DrugCLIP '+DATA.g2, dcB, 16); trackRefs['dcB'] = dcB;

  // Register track rows into toggle groups
  trackGroupRows = {
    structure: [
      alndARow.row, alndBRow.row,
      tedARow.row, tedBRow.row,
      disorderARow.row, disorderBRow.row,
      helixARow.row, helixBRow.row,
      strandARow.row, strandBRow.row,
    ],
    conservation: [
      plmaARow.row, ...plmaSubRowsA,
      plmaBRow.row, ...plmaSubRowsB,
    ],
    druggability: [
      cavARow.row, ...cavSubRowsA,
      cavBRow.row, ...cavSubRowsB,
      dcARow.row, dcBRow.row,
    ],
  };
  // Apply current group visibility state
  Object.entries(trackGroupRows).forEach(([group, rows]) => {
    const visible = trackGroupState[group];
    rows.forEach(r => { r.style.display = visible ? '' : 'none'; });
  });

  requestAnimationFrame(()=>{
    const allTracks = [
      nav, hlTrackA, hlTrackB, seqA, seqB,
      domA, domB, tedA, tedB,
      plmaCatA, plmaSharedA, plmaPairA, plmaSpecA, plmaFamA,
      plmaCatB, plmaSharedB, plmaPairB, plmaSpecB, plmaFamB,
      amA, dam, amB,
      alndA, alndB,
      disorderA, disorderB, helixA, helixB, strandA, strandB,
      cavA, cavStrongA, cavMediumA, cavWeakA,
      cavB, cavStrongB, cavMediumB, cavWeakB,
      dcA, dcB
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

    [amA, amB, dam, cavA, cavB, plmaCatA, plmaCatB, hlTrackA, hlTrackB].forEach(track => {
      track.setAttribute('shape','rectangle');
    });
    hlTrackA.data = [];
    hlTrackB.data = [];
    [plmaSharedA, plmaPairA, plmaSpecA, plmaFamA, plmaSharedB, plmaPairB, plmaSpecB, plmaFamB].forEach(track => {
      track.setAttribute('shape','rectangle');
    });
    helixA.setAttribute('shape','helix'); helixB.setAttribute('shape','helix');
    strandA.setAttribute('shape','strand'); strandB.setAttribute('shape','strand');

    [domA, disorderA, tedA, cavStrongA, cavMediumA, cavWeakA, dcA, domB, disorderB, tedB, cavStrongB, cavMediumB, cavWeakB, dcB].forEach(track => {
      track.setAttribute('shape','roundRectangle');
      track.setAttribute('show-label','');
    });

    const allSsA = sanitizeRects(DATA.ssA_alnRects||[], alnLen);
    helixA.data = allSsA.filter(r => r.ss_type === 'helix');
    strandA.data = allSsA.filter(r => r.ss_type === 'strand');
    const allSsB = sanitizeRects(DATA.ssB_alnRects||[], alnLen);
    helixB.data = allSsB.filter(r => r.ss_type === 'helix');
    strandB.data = allSsB.filter(r => r.ss_type === 'strand');

    domA.data       = sanitizeRects(DATA.domA_alnRects||[], alnLen);          domA._originalData = [...domA.data];
    disorderA.data  = sanitizeRects(DATA.disorderA_alnRects||[], alnLen);     disorderA._originalData = [...disorderA.data];
    tedA.data       = sanitizeRects(DATA.tedA_alnRects||[], alnLen);          tedA._originalData = [...tedA.data];
    const allCavA = sanitizeRects(DATA.cavA_alnRects||[], alnLen);
    cavA._originalData = [...allCavA];
    cavStrongA._originalData = allCavA.filter(r => (r.druggability||'').toLowerCase() === 'strong');
    cavMediumA._originalData = allCavA.filter(r => (r.druggability||'').toLowerCase() === 'medium');
    cavWeakA._originalData = allCavA.filter(r => (r.druggability||'').toLowerCase() === 'weak');
    dcA.data        = sanitizeRects(DATA.dcA_alnRects||[], alnLen);           dcA._originalData = [...dcA.data];

    domB.data       = sanitizeRects(DATA.domB_alnRects||[], alnLen);          domB._originalData = [...domB.data];
    disorderB.data  = sanitizeRects(DATA.disorderB_alnRects||[], alnLen);     disorderB._originalData = [...disorderB.data];
    tedB.data       = sanitizeRects(DATA.tedB_alnRects||[], alnLen);          tedB._originalData = [...tedB.data];
    const allCavB = sanitizeRects(DATA.cavB_alnRects||[], alnLen);
    cavB._originalData = [...allCavB];
    cavStrongB._originalData = allCavB.filter(r => (r.druggability||'').toLowerCase() === 'strong');
    cavMediumB._originalData = allCavB.filter(r => (r.druggability||'').toLowerCase() === 'medium');
    cavWeakB._originalData = allCavB.filter(r => (r.druggability||'').toLowerCase() === 'weak');
    dcB.data        = sanitizeRects(DATA.dcB_alnRects||[], alnLen);           dcB._originalData = [...dcB.data];

    // PLMA category tracks: overview (priority heatmap) + sub-tracks per category
    plmaCatA.data = buildPlmaCategoryRects('A');
    plmaSharedA.data = buildPlmaCategoryRects('A', 'shared_with_family');
    plmaPairA.data = buildPlmaCategoryRects('A', 'pair_exclusive');
    plmaSpecA.data = buildPlmaCategoryRects('A', 'specific_a');
    plmaFamA.data = buildPlmaCategoryRects('A', 'a_with_family');

    plmaCatB.data = buildPlmaCategoryRects('B');
    plmaSharedB.data = buildPlmaCategoryRects('B', 'shared_with_family');
    plmaPairB.data = buildPlmaCategoryRects('B', 'pair_exclusive');
    plmaSpecB.data = buildPlmaCategoryRects('B', 'specific_b');
    plmaFamB.data = buildPlmaCategoryRects('B', 'b_with_family');

    // Aligned substitution type tracks
    alndA.setAttribute('shape', 'rectangle');
    alndB.setAttribute('shape', 'rectangle');
    alndA.data = sanitizeRects(buildAlignedRects('A'), alnLen);
    alndB.data = sanitizeRects(buildAlignedRects('B'), alnLen);

    applyAmMode(amMode);
    applyCavityFilter(); // Apply default druggability filter (medium+)
  });

  function attachDomainClick(track, chain) {
    const domMap = (chain === chainIdA) ? domByUidA : domByUidB;
    // Click only — find feature at clicked alignment position
    track.addEventListener('click', (ev)=>{
      const trackData = track.data || [];
      if (!trackData.length) return;
      const mgrEl = document.getElementById('mgr');
      const dispStart = parseInt(mgrEl?.getAttribute('display-start') || track.getAttribute('display-start') || '1');
      const dispEnd = parseInt(mgrEl?.getAttribute('display-end') || track.getAttribute('display-end') || String(alnLen));
      const rect = track.getBoundingClientRect();
      const fraction = ev.clientX - rect.left;
      const trackWidth = rect.width;
      if (trackWidth <= 0) return;
      const pos = Math.round(dispStart + (fraction / trackWidth) * (dispEnd - dispStart));
      for (const f of trackData) {
        if (!f.id) continue;
        const s = f.start || f.x || f.begin || 1;
        const e = f.end || f.to || s;
        if (pos >= s && pos <= e) {
          const dom = domMap[f.id];
          if (dom) { toggleFeature(dom, chain); return; }
        }
      }
    });
  }

  attachDomainClick(domA, chainIdA);
  attachDomainClick(disorderA, chainIdA);
  attachDomainClick(tedA, chainIdA);
  attachDomainClick(cavA, chainIdA);
  attachDomainClick(cavStrongA, chainIdA);
  attachDomainClick(cavMediumA, chainIdA);
  attachDomainClick(cavWeakA, chainIdA);
  attachDomainClick(dcA, chainIdA);
  attachDomainClick(domB, chainIdB);
  attachDomainClick(disorderB, chainIdB);
  attachDomainClick(tedB, chainIdB);
  attachDomainClick(cavB, chainIdB);
  attachDomainClick(cavStrongB, chainIdB);
  attachDomainClick(cavMediumB, chainIdB);
  attachDomainClick(cavWeakB, chainIdB);
  attachDomainClick(dcB, chainIdB);
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
  // Regions & Druggability
  setupCollapsibleSection('regionsCollapseBtn', 'regionsBody', 'regionsSection');
  // Family group wrapper — with re-init constellation on first expand
  setupCollapsibleSection('familyGroupCollapseBtn', 'familyGroupBody', 'familyGroup');
  let familyConstellationNeedsInit = true;
  const familyGroupBtn = document.getElementById('familyGroupCollapseBtn');
  if (familyGroupBtn) {
    familyGroupBtn.addEventListener('click', () => {
      const body = document.getElementById('familyGroupBody');
      if (body && !body.classList.contains('collapsed') && familyConstellationNeedsInit) {
        familyConstellationNeedsInit = false;
        // Delay to allow layout reflow after un-collapsing
        requestAnimationFrame(() => { initFamilyConstellation(); });
      }
    }, { passive: true });
  }
}

// Apply default collapse states: expand key sections, collapse secondary ones
function applyDefaultCollapseStates() {
  const collapseByDefault = [
    { btn: 'pdbeCollapseBtn', body: 'pdbeCardBody', section: 'pdbeCard' },
  ];
  collapseByDefault.forEach(({ btn, body, section }) => {
    const b = document.getElementById(body);
    const bt = document.getElementById(btn);
    const s = section ? document.getElementById(section) : null;
    if (b && !b.classList.contains('collapsed')) {
      b.classList.add('collapsed');
      if (bt) { bt.setAttribute('aria-expanded', 'false'); bt.textContent = 'Expand'; }
      if (s) s.classList.add('is-collapsed');
    }
  });
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

  // Regions section - hide if no domains at all
  const hasDomainsA = DATA && DATA.domainsA && DATA.domainsA.length > 0;
  const hasDomainsB = DATA && DATA.domainsB && DATA.domainsB.length > 0;
  if (!hasDomainsA && !hasDomainsB) {
    hideSection('regionsSection');
  } else {
    showSection('regionsSection');
  }

  // Family features section - hide if no family features data or family of 2
  const hasFamilyFeatures = SUMMARY && SUMMARY.family_features && Object.keys(SUMMARY.family_features).length > 0;
  if (!hasFamilyFeatures || isFamily2) {
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

  // Family group - hide if all inner sections are hidden
  const familyNavHidden = document.getElementById('familyNav')?.classList.contains('section-hidden') || document.getElementById('familyNav')?.style.display === 'none';
  const familyFeatHidden = document.getElementById('familyFeaturesSection')?.classList.contains('section-hidden');
  const simSearchHidden = document.getElementById('similaritySearchSection')?.classList.contains('section-hidden');
  if (familyNavHidden && familyFeatHidden && simSearchHidden) {
    hideSection('familyGroup');
  } else {
    showSection('familyGroup');
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
  if (druggabilityFilter === 'strong') return dg === 'strong';
  if (druggabilityFilter === 'medium+') return dg === 'strong' || dg === 'medium';

  return true;
}

function shouldShowCavityRect(rect) {
  // Filter cavity rectangles based on druggability
  if (!rect.druggability) return true; // No druggability info, show it

  const dg = (rect.druggability || '').toLowerCase();

  if (druggabilityFilter === 'all') return true;
  if (druggabilityFilter === 'strong') return dg === 'strong';
  if (druggabilityFilter === 'medium+') return dg === 'strong' || dg === 'medium';

  return true;
}

function buildCavityOverview(rects, alnLen) {
  // Build per-residue heatmap: each position gets the strongest druggability color
  const PRI = { 'strong': 2, 'medium': 1, 'weak': 0 };
  const COLORS = { 2: '#e65100', 1: '#ff9800', 0: '#ffcc80' };

  const best = new Int8Array(alnLen + 1).fill(-1);
  for (const r of rects) {
    const s = r.x || r.start || 1;
    const e = r.end || r.to || s;
    const pri = PRI[(r.druggability||'').toLowerCase()] ?? -1;
    if (pri < 0) continue;
    for (let p = s; p <= e && p <= alnLen; p++) {
      if (pri > best[p]) best[p] = pri;
    }
  }

  // Merge adjacent positions with same priority into runs
  const out = [];
  let runStart = -1, runPri = -1;
  for (let p = 1; p <= alnLen + 1; p++) {
    const cur = p <= alnLen ? best[p] : -1;
    if (cur === runPri && cur >= 0) continue;
    if (runStart > 0 && runPri >= 0) {
      out.push({ x: runStart, start: runStart, begin: runStart, end: p - 1, to: p - 1, color: COLORS[runPri], opacity: 0.8 });
    }
    runStart = p; runPri = cur;
  }
  return out;
}

// PLMA category colors: vivid for key categories, muted for background
const PLMA_CAT_COLORS = {
  specific_a:         '#EF5350',  // red — gene A only
  specific_b:         '#AB47BC',  // purple — gene B only
  pair_exclusive:     '#26A69A',  // teal — both paralogs only
  a_with_family:      '#FF7043',  // deep orange — gene A + family
  b_with_family:      '#7E57C2',  // deep purple — gene B + family
  shared_with_family: '#78909C',  // blue-grey — both + family
  family_only:        '#BDBDBD',  // gray — other family only
};
// Priority for overview: gene-specific > gene+family > pair > shared > family-only
const PLMA_CAT_PRI = { specific_a:5, specific_b:5, a_with_family:4, b_with_family:4, pair_exclusive:3, shared_with_family:2, family_only:1 };

// Build per-residue PLMA category rects for a gene, mapped to pairwise alignment columns.
// If filterCat is given, only include blocks of that category (for sub-tracks).
function buildPlmaCategoryRects(side, filterCat) {
  if (!PLMA_DATA || !DATA) return [];
  const aln = side === 'A' ? (DATA.qaln || '') : (DATA.taln || '');
  const seqNum = side === 'A' ? PLMA_DATA.gene_a_seq : PLMA_DATA.gene_b_seq;
  if (!seqNum || !aln) return [];

  const posToCol = {};
  let pos = 0;
  for (let col = 0; col < aln.length; col++) {
    if (aln[col] !== '-') { pos++; posToCol[pos] = col + 1; }
  }
  const alnLen = aln.length;

  const colCat = new Array(alnLen + 1).fill(null);
  const colPri = new Int8Array(alnLen + 1).fill(-1);

  for (const block of (PLMA_DATA.blocks || [])) {
    const p = block.positions?.[seqNum];
    if (!p) continue;
    const cat = block.category;
    if (filterCat && cat !== filterCat) continue;
    const pri = filterCat ? 1 : (PLMA_CAT_PRI[cat] ?? -1);
    for (let r = p.start; r <= p.end; r++) {
      const c = posToCol[r];
      if (!c || c > alnLen) continue;
      if (pri > colPri[c]) { colPri[c] = pri; colCat[c] = cat; }
    }
  }

  const out = [];
  let runStart = -1, runCat = null;
  for (let c = 1; c <= alnLen + 1; c++) {
    const cur = c <= alnLen ? colCat[c] : null;
    if (cur === runCat && cur !== null) continue;
    if (runStart > 0 && runCat) {
      out.push({ x:runStart, start:runStart, begin:runStart, end:c-1, to:c-1,
                 color: PLMA_CAT_COLORS[runCat] || '#ccc', opacity:0.85 });
    }
    runStart = c; runCat = cur;
  }
  return out;
}

function applyCavityFilter() {
  const alnLen = Math.max(1, (DATA.qaln||'').length);

  ['A', 'B'].forEach(side => {
    const ov = trackRefs['cav' + side];
    const strong = trackRefs['cavStrong' + side];
    const medium = trackRefs['cavMedium' + side];
    const weak = trackRefs['cavWeak' + side];

    if (ov && ov._originalData) {
      const filtered = ov._originalData.filter(shouldShowCavityRect);
      ov.data = buildCavityOverview(filtered, alnLen);
    }
    if (strong && strong._originalData) {
      strong.data = strong._originalData.filter(shouldShowCavityRect);
    }
    if (medium && medium._originalData) {
      medium.data = medium._originalData.filter(shouldShowCavityRect);
    }
    if (weak && weak._originalData) {
      weak.data = weak._originalData.filter(shouldShowCavityRect);
    }
  });
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
  DATA.qposToCol = qposToCol;
  DATA.tposToCol = tposToCol;

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

  assignDomainColors(DATA.domainsA, DATA.domainsB);

  // Recalculate domain/disorder/cavity/drugclip rectangles
  if (DATA.domainsA) {
    const doms = DATA.domainsA.filter(d => d.type !== 'CAV' && d.type !== 'Cavity' && d.type !== 'DrugCLIP');
    const cavs = DATA.domainsA.filter(d => d.type === 'CAV' || d.type === 'Cavity');
    const dcs  = DATA.domainsA.filter(d => d.type === 'DrugCLIP');

    DATA.domA_alnRects = remapFeaturesToAlignment(doms, qposToCol, '#2ca02c');
    DATA.cavA_alnRects = remapFeaturesToAlignment(cavs, qposToCol, '#ff7d45');
    DATA.dcA_alnRects  = remapFeaturesToAlignment(dcs,  qposToCol, '#c62828');
  }

  if (DATA.domainsB) {
    const doms = DATA.domainsB.filter(d => d.type !== 'CAV' && d.type !== 'Cavity' && d.type !== 'DrugCLIP');
    const cavs = DATA.domainsB.filter(d => d.type === 'CAV' || d.type === 'Cavity');
    const dcs  = DATA.domainsB.filter(d => d.type === 'DrugCLIP');

    DATA.domB_alnRects = remapFeaturesToAlignment(doms, tposToCol, '#2ca02c');
    DATA.cavB_alnRects = remapFeaturesToAlignment(cavs, tposToCol, '#ff7d45');
    DATA.dcB_alnRects  = remapFeaturesToAlignment(dcs,  tposToCol, '#c62828');
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
  // Gene panel headers
  const hA = document.getElementById('regionsPanelAHeader');
  const hB = document.getElementById('regionsPanelBHeader');
  if (hA) hA.textContent = DATA.g1 || 'Gene A';
  if (hB) hB.textContent = DATA.g2 || 'Gene B';

  // Split domains into structural regions vs druggability
  const isStructural = d => d.type !== 'Cavity' && d.raw_type !== 'Cavity' && d.type !== 'DrugCLIP' && d.raw_type !== 'DrugCLIP';
  const isDrug = d => !isStructural(d);

  const regA = document.getElementById('regionsTableA'); if (regA) regA.innerHTML = '';
  const regB = document.getElementById('regionsTableB'); if (regB) regB.innerHTML = '';
  const drugA = document.getElementById('drugTableA'); if (drugA) drugA.innerHTML = '';
  const drugB = document.getElementById('drugTableB'); if (drugB) drugB.innerHTML = '';

  // OpenChemLib availability
  const _oclOk = typeof OCL !== 'undefined' && typeof OCL.Molecule !== 'undefined';
  function renderSmilesToSvg(container, smiles) {
    if (!_oclOk || !smiles || !container) return;
    try {
      const mol = OCL.Molecule.fromSmiles(smiles);
      const svgStr = mol.toSVG(200, 150);
      const wrap = document.createElement('div');
      wrap.style.cssText = 'width:200px;height:150px;border:1px solid #eee;border-radius:6px;background:#fff;flex-shrink:0;overflow:hidden';
      wrap.innerHTML = svgStr;
      container.prepend(wrap);
    } catch(e) {
      const fallback = document.createElement('div');
      fallback.style.cssText = 'width:200px;height:150px;border:1px solid #eee;border-radius:6px;background:#fff;display:flex;align-items:center;justify-content:center;font-size:9px;color:#999;text-align:center;padding:8px;word-break:break-all;flex-shrink:0';
      fallback.textContent = smiles.length > 80 ? smiles.slice(0, 80) + '...' : smiles;
      container.prepend(fallback);
    }
  }
  function drugDbLink(oid) {
    if (!oid) return null;
    if (oid.startsWith('ZINC')) return { url: `https://zinc.docking.org/substances/${oid}/`, label: 'ZINC' };
    return null;
  }

  // Structural region row
  const addRegionRow = (tb, d, chain) => {
    const tr = document.createElement('tr'); tr.className = 'clickable';
    tr.setAttribute('data-uid', d.uid || ''); tr.setAttribute('data-chain', chain);
    const tdSel = document.createElement('td');
    const cb = document.createElement('input'); cb.type = 'checkbox';
    cb.addEventListener('change', ev => { ev.stopPropagation(); toggleFeature(d, chain); }, {passive:true});
    tdSel.append(cb);
    const tdName = document.createElement('td');
    tdName.textContent = d.label || d.name || d.type || '';
    if (d.color) tdName.style.cssText = `border-left:3px solid ${d.color};padding-left:6px`;
    const tdRange = document.createElement('td'); tdRange.textContent = `${d.start}-${d.end}`;
    tr.append(tdSel, tdName, tdRange);
    tr.addEventListener('click', ev => {
      if (ev.target.closest('input[type="checkbox"]')) return;
      toggleFeature(d, chain);
    }, {passive:true});
    tb.append(tr);
  };

  // Druggability row (cavity or DrugCLIP)
  const addDrugRow = (tb, d, chain) => {
    const tr = document.createElement('tr'); tr.className = 'clickable';
    tr.setAttribute('data-uid', d.uid || ''); tr.setAttribute('data-chain', chain);
    const tdSel = document.createElement('td');
    const cb = document.createElement('input'); cb.type = 'checkbox';
    cb.addEventListener('change', ev => { ev.stopPropagation(); toggleFeature(d, chain); }, {passive:true});
    tdSel.append(cb);
    const tdName = document.createElement('td');
    tdName.textContent = d.label || d.name || d.type || '';
    const tdRange = document.createElement('td'); tdRange.textContent = `${d.start}-${d.end}`;
    const tdScore = document.createElement('td');

    if (d.type === 'DrugCLIP' || d.raw_type === 'DrugCLIP') {
      const nTotal = d.n_total_hits || 0;
      const drugs = d.top_drugs || [];
      if (drugs.length > 0) {
        const badge = document.createElement('span');
        badge.style.cssText = 'background:#c62828;color:#fff;padding:1px 8px;border-radius:10px;font-size:10px;cursor:pointer;white-space:nowrap';
        badge.textContent = `${nTotal} hits ▾`;
        tdScore.append(badge);
        // Expandable detail row below
        const expandTr = document.createElement('tr'); expandTr.className = 'dc-expand-row';
        expandTr.style.display = 'none';
        const expandTd = document.createElement('td'); expandTd.colSpan = 4;
        expandTd.style.cssText = 'padding:8px 10px';
        expandTr.append(expandTd);

        const toggleExpand = (ev) => {
          ev.stopPropagation();
          const open = expandTr.style.display !== 'none';
          expandTr.style.display = open ? 'none' : '';
          badge.textContent = open ? `${nTotal} hits ▾` : `${nTotal} hits ▴`;
          if (!expandTd._built) {
            expandTd._built = true;
            buildDrugExpansion(expandTd, drugs, d, chain, _oclOk, renderSmilesToSvg, drugDbLink);
          }
        };
        badge.addEventListener('click', toggleExpand);
        tr.addEventListener('click', ev => {
          if (ev.target.closest('input[type="checkbox"]') || ev.target.closest('span')) return;
          toggleFeature(d, chain);
        }, {passive:true});
        tr.append(tdSel, tdName, tdRange, tdScore);
        tb.append(tr, expandTr);
        return;
      } else {
        tdScore.textContent = 'No hits';
      }
    } else if (d.type === 'Cavity' || d.raw_type === 'Cavity') {
      const dg = d.druggability || '';
      const ds = d.drug_score || d.drugscore || '';
      tdScore.textContent = (ds || dg) ? `${ds || ''} ${dg || ''}`.trim() : '';
      if (dg === 'strong') tdScore.style.color = '#16a34a';
      else if (dg === 'medium') tdScore.style.color = '#d97706';
    }

    tr.append(tdSel, tdName, tdRange, tdScore);
    tr.addEventListener('click', ev => {
      if (ev.target.closest('input[type="checkbox"]')) return;
      toggleFeature(d, chain);
    }, {passive:true});
    tb.append(tr);
  };

  // Fill structural regions
  (DATA.domainsA||[]).filter(d => isStructural(d) && shouldShowDomain(d)).forEach(d => addRegionRow(regA, d, chainIdA));
  (DATA.domainsB||[]).filter(d => isStructural(d) && shouldShowDomain(d)).forEach(d => addRegionRow(regB, d, chainIdB));
  // Fill druggability (cavities + DrugCLIP)
  (DATA.domainsA||[]).filter(d => isDrug(d) && shouldShowDomain(d)).forEach(d => addDrugRow(drugA, d, chainIdA));
  (DATA.domainsB||[]).filter(d => isDrug(d) && shouldShowDomain(d)).forEach(d => addDrugRow(drugB, d, chainIdB));

  // Hide DxD section if no domain pairs
  const dxd = document.getElementById('dxdSection');
  if (dxd) dxd.style.display = (Array.isArray(DATA.domPairs) && DATA.domPairs.length) ? '' : 'none';

  // Hide regions section if no data at all
  const hasAny = (DATA.domainsA||[]).length > 0 || (DATA.domainsB||[]).length > 0;
  const regSection = document.getElementById('regionsSection');
  if (regSection) regSection.classList.toggle('section-hidden', !hasAny);

  renderTableSelections();
}

// Build the expandable drug hits content for a DrugCLIP pocket
function buildDrugExpansion(container, drugs, pocket, chain, _oclOk, renderSmilesToSvg, drugDbLink) {
  const tbl = document.createElement('table');
  tbl.style.cssText = 'width:100%;border-collapse:collapse;font-size:11px';
  tbl.innerHTML = '<thead><tr style="border-bottom:1px solid #ddd;background:#f9f7f2"><th style="text-align:left;padding:3px 6px">#</th><th style="text-align:left;padding:3px 6px">SMILES</th><th style="text-align:left;padding:3px 6px">ID</th><th style="text-align:right;padding:3px 6px">Score</th></tr></thead>';
  const tbody = document.createElement('tbody');
  drugs.forEach((drug, idx) => {
    const tr = document.createElement('tr');
    tr.style.cssText = 'border-bottom:1px solid #f0f0f0;cursor:pointer';
    tr.addEventListener('mouseenter', () => { if (!tr._expanded) tr.style.background = '#fff5f5'; });
    tr.addEventListener('mouseleave', () => { if (!tr._expanded) tr.style.background = ''; });
    const smi = drug.smiles || '';
    const smiShort = smi.length > 30 ? smi.slice(0, 30) + '...' : smi;
    const oid = drug.oid || '';
    const oidShort = oid.length > 18 ? oid.slice(0, 18) + '...' : oid;
    const dbLnk = drugDbLink(oid);
    const oidCell = dbLnk
      ? `<a href="${dbLnk.url}" target="_blank" style="color:#c62828" title="${oid}" onclick="event.stopPropagation()">${oidShort}</a>`
      : `<span title="${oid}">${oidShort}</span>`;
    tr.innerHTML = `<td style="padding:3px 6px">${idx + 1}</td><td style="padding:3px 6px;font-family:monospace;font-size:10px" title="${smi}">${smiShort}</td><td style="padding:3px 6px;font-size:10px">${oidCell}</td><td style="padding:3px 6px;text-align:right;font-weight:600">${drug.score ? drug.score.toFixed(2) : '–'}</td>`;

    const detailTr = document.createElement('tr');
    detailTr.style.display = 'none';
    const detailTd = document.createElement('td'); detailTd.colSpan = 4;
    detailTd.style.cssText = 'padding:8px 6px;background:#fafafa;border-bottom:1px solid #e0e0e0';
    detailTr.append(detailTd);

    tr.addEventListener('click', () => {
      const wasExpanded = tr._expanded;
      tbody.querySelectorAll('tr[data-expanded="true"]').forEach(r => {
        r._expanded = false; r.removeAttribute('data-expanded'); r.style.background = '';
        const next = r.nextElementSibling;
        if (next && next.querySelector('td[colspan]')) next.style.display = 'none';
      });
      if (!wasExpanded) {
        tr._expanded = true; tr.setAttribute('data-expanded', 'true'); tr.style.background = '#fff0f0';
        detailTr.style.display = '';
        if (!detailTd._built) {
          detailTd._built = true;
          const wrap = document.createElement('div');
          wrap.style.cssText = 'display:flex;gap:12px;align-items:flex-start';
          if (_oclOk && smi) renderSmilesToSvg(wrap, smi);
          const info = document.createElement('div'); info.style.cssText = 'font-size:11px;min-width:0';
          const dbLink = drugDbLink(oid);
          const idHtml = dbLink
            ? `<a href="${dbLink.url}" target="_blank" style="color:#c62828;text-decoration:underline">${oid}</a> <span style="color:#999">(${dbLink.label})</span>`
            : oid;
          info.innerHTML = `<div style="margin-bottom:4px"><strong>Full SMILES:</strong></div><div style="font-family:monospace;font-size:10px;word-break:break-all;background:#fff;padding:4px 6px;border:1px solid #eee;border-radius:4px;margin-bottom:6px">${smi}</div><div><strong>ID:</strong> ${idHtml}</div><div><strong>Score:</strong> ${drug.score ? drug.score.toFixed(4) : '–'}</div>`;
          const hlBtn = document.createElement('button');
          hlBtn.style.cssText = 'margin-top:6px;background:#c62828;color:#fff;border:none;border-radius:4px;padding:3px 10px;cursor:pointer;font-size:11px';
          hlBtn.textContent = 'Highlight pocket in 3D';
          hlBtn.addEventListener('click', ev => { ev.stopPropagation(); toggleFeature(pocket, chain); });
          info.append(hlBtn);
          wrap.append(info);
          detailTd.append(wrap);
        }
      }
    });
    tbody.append(tr, detailTr);
  });
  tbl.append(tbody);
  container.append(tbl);
}

async function fillDomPairs(){
  const tb = document.querySelector('#domPairs tbody'); tb.innerHTML='';
  if (!Array.isArray(DATA.domPairs) || !DATA.domPairs.length){
    tb.innerHTML = '<tr><td colspan="5" class="small">No domain sub-alignments computed.</td></tr>';
    return;
  }
  DATA.domPairs.forEach((r, idx)=>{
    const tr=document.createElement('tr'); tr.className='clickable';
    tr.innerHTML = `<td>${r.Aname} ${r.Arng}</td><td>${r.Bname} ${r.Brng}</td><td>${r.fident!=null?r.fident.toFixed(1)+'%':'–'}</td><td>${r.tm!=null?r.tm.toFixed(3):'–'}</td><td>${r.damPct!=null?r.damPct.toFixed(1)+'%':'–'}</td>`;
    tr.addEventListener('click', async ()=>{
      await reloadViewerWith(r.pdb64);
      setTmScoreDisplay(r.tm);
      const title = `Domain: ${r.Aname} ${r.Arng} × ${r.Bname} ${r.Brng}`;
      document.getElementById('contextTitle').textContent = title;

      selection.clear();
      const domA = (DATA.domainsA||[]).find(d => (d.label===r.Aname) || (`${d.start}-${d.end}`===r.Arng));
      const domB = (DATA.domainsB||[]).find(d => (d.label===r.Bname) || (`${d.start}-${d.end}`===r.Brng));
      if (domA) selection.set(selectionKey(chainIdA, domA.uid), {
        id: domA.uid, chain: chainIdA,
        start: parseInt(domA.start,10), end: parseInt(domA.end,10),
        color: domA.color||'#ffdb13', name: domA.label||domA.name
      });
      if (domB) selection.set(selectionKey(chainIdB, domB.uid), {
        id: domB.uid, chain: chainIdB,
        start: parseInt(domB.start,10), end: parseInt(domB.end,10),
        color: domB.color||'#ff7d45', name: domB.label||domB.name
      });
      await renderSelections();
      syncHighlightDropdowns();
      const alnSel = document.getElementById('alnSelector');
      if (alnSel) alnSel.value = String(idx);
    }, {passive:true});
    tb.append(tr);
  });
}

// fillDrugHits is now absorbed into fillDomainTables
function fillDrugHits() { /* no-op — drug hits are rendered inline in fillDomainTables */ }

let currentColorMode = 'uniform';

// --- Client-side PDB B-factor modification for Molstar coloring ---
const _pdbBfCache = {};

function modifyPdbBfactors(bfMapA, bfMapB, defaultBf) {
  if (!PDB64_FULL) return null;
  const pdbText = atob(PDB64_FULL);
  const lines = pdbText.split('\n');
  const out = [];
  for (const line of lines) {
    if (line.length >= 66 && (line.startsWith('ATOM') || line.startsWith('HETATM'))) {
      const chain = line[21];
      const resSeq = parseInt(line.substring(22, 26).trim(), 10);
      const map = chain === 'A' ? bfMapA : bfMapB;
      const bf = (map && resSeq in map) ? map[resSeq] : defaultBf;
      out.push(line.substring(0, 60) + bf.toFixed(2).padStart(6) + line.substring(66));
    } else {
      out.push(line);
    }
  }
  return btoa(out.join('\n'));
}

function alnColToResMap(aln) {
  const m = {};
  let pos = 0;
  for (let col = 0; col < aln.length; col++) {
    if (aln[col] !== '-') { pos++; m[col + 1] = pos; }
  }
  return m;
}

function buildPlddtBfactorMaps() {
  if (!DATA) return null;
  const mapA = {}, mapB = {};
  let plddtA = DATA.plddtA || [], plddtB = DATA.plddtB || [];
  // Fallback: extract pLDDT from AF PDB B-factors when DATA arrays are empty
  if (!plddtA.length && window.PDB64_A) plddtA = extractPlddtFromPdb(window.PDB64_A);
  if (!plddtB.length && window.PDB64_B) plddtB = extractPlddtFromPdb(window.PDB64_B);
  if (!plddtA.length && !plddtB.length) return null;
  // Bin to match AlphaFold official cutoffs: <=50→0, 51-70→1, 71-90→2, >90→3
  const bin = v => v <= 50 ? 0 : v <= 70 ? 1 : v <= 90 ? 2 : 3;
  plddtA.forEach((v, i) => { if (typeof v === 'number') mapA[i + 1] = bin(v); });
  plddtB.forEach((v, i) => { if (typeof v === 'number') mapB[i + 1] = bin(v); });
  return { A: mapA, B: mapB };
}

function buildAmBfactorMaps() {
  if (!DATA) return null;
  const mapA = {}, mapB = {};
  (DATA.bfactorsA || []).forEach((v, i) => { if (typeof v === 'number') mapA[i + 1] = v * 100; });
  (DATA.bfactorsB || []).forEach((v, i) => { if (typeof v === 'number') mapB[i + 1] = v * 100; });
  return { A: mapA, B: mapB };
}

function buildDamBfactorMaps() {
  if (!DATA) return null;
  const mapA = {}, mapB = {};
  const qaln = DATA.qaln || '', taln = DATA.taln || '';
  const colToA = alnColToResMap(qaln), colToB = alnColToResMap(taln);
  const bfA = DATA.bfactorsA || [], bfB = DATA.bfactorsB || [];
  for (let col = 1; col <= qaln.length; col++) {
    const pa = colToA[col], pb = colToB[col];
    if (pa && pb && pa <= bfA.length && pb <= bfB.length) {
      const va = bfA[pa - 1], vb = bfB[pb - 1];
      if (typeof va === 'number' && typeof vb === 'number') {
        const delta = Math.abs(va - vb) * 100;
        mapA[pa] = delta;
        mapB[pb] = delta;
      }
    }
  }
  return { A: mapA, B: mapB };
}

// Amino acid biochemical class groups for substitution classification
const AA_CLASSES = {G:0,A:0,V:0,L:0,I:0,P:0,M:0,F:1,W:1,Y:1,S:2,T:2,C:2,N:2,Q:2,K:3,R:3,H:3,D:4,E:4};
function aaSubType(a, b) {
  // Returns 0=gap, 0.25=identical, 0.5=conservative (same class), 1=radical (different class)
  if (!a || a === '-' || !b || b === '-') return 0;
  if (a.toUpperCase() === b.toUpperCase()) return 0.25; // identical → grey
  const ca = AA_CLASSES[a.toUpperCase()], cb = AA_CLASSES[b.toUpperCase()];
  if (ca === undefined || cb === undefined) return 0.5;
  return ca === cb ? 0.5 : 1; // same class=conservative(teal), different=radical(green)
}

function buildAlignedBfactorMaps() {
  if (!DATA) return null;
  const mapA = {}, mapB = {};
  const qaln = DATA.qaln || '', taln = DATA.taln || '';
  const colToA = alnColToResMap(qaln), colToB = alnColToResMap(taln);
  for (let col = 1; col <= qaln.length; col++) {
    const pa = colToA[col], pb = colToB[col];
    const qa = qaln[col - 1], ta = taln[col - 1];
    const subType = (pa && pb) ? aaSubType(qa, ta) : 0;
    if (pa) mapA[pa] = subType;
    if (pb) mapB[pb] = subType;
  }
  return { A: mapA, B: mapB };
}

// Build alignment coloring rects for Nightingale tracks
// Returns array of {start, end, color} in 1-indexed alignment column space
function buildAlignedRects(chain) {
  if (!DATA) return [];
  const qaln = DATA.qaln || '', taln = DATA.taln || '';
  const myAln = chain === 'A' ? qaln : taln;
  const otherAln = chain === 'A' ? taln : qaln;
  const SUB_COLORS = {0: '#8B1A1A', 0.25: '#AAAAAA', 0.5: '#009688', 1: '#4CAF50'};
  const rects = [];
  let curType = null, curStart = null;
  for (let col = 0; col < myAln.length; col++) {
    const myAA = myAln[col];
    if (myAA === '-') {
      if (curStart !== null) {
        rects.push({start: curStart, end: col, color: SUB_COLORS[curType]});
        curStart = null; curType = null;
      }
      continue;
    }
    const otherAA = otherAln[col];
    const subType = (otherAA && otherAA !== '-') ? aaSubType(myAA, otherAA) : 0;
    if (subType !== curType) {
      if (curStart !== null) rects.push({start: curStart, end: col, color: SUB_COLORS[curType]});
      curStart = col + 1; curType = subType;
    }
  }
  if (curStart !== null) rects.push({start: curStart, end: myAln.length, color: SUB_COLORS[curType]});
  return rects;
}

function buildDomainBfactorMaps() {
  if (!DATA) return null;
  // Assign a unique index per domain name (consistent across A and B)
  const nameToIdx = {};
  let idx = 1;
  const allDoms = [...(DATA.domainsA||[]), ...(DATA.domainsB||[])];
  for (const d of allDoms) {
    if (d.type !== 'Domain' && d.raw_type !== 'DOMAIN') continue;
    const name = d.label || d.name || d.type || 'unknown';
    if (!(name in nameToIdx)) nameToIdx[name] = idx++;
  }
  // Store for theme generation
  window._domainColorNames = nameToIdx;
  const mapA = {}, mapB = {};
  for (const d of (DATA.domainsA || [])) {
    if (d.type === 'Domain' || d.raw_type === 'DOMAIN') {
      const name = d.label || d.name || d.type || 'unknown';
      const val = nameToIdx[name] || 1;
      for (let r = d.start; r <= d.end; r++) mapA[r] = val;
    }
  }
  for (const d of (DATA.domainsB || [])) {
    if (d.type === 'Domain' || d.raw_type === 'DOMAIN') {
      const name = d.label || d.name || d.type || 'unknown';
      const val = nameToIdx[name] || 1;
      for (let r = d.start; r <= d.end; r++) mapB[r] = val;
    }
  }
  return { A: mapA, B: mapB };
}

function rectsToResMap(rects, colToRes, valueFn) {
  const map = {};
  for (const r of rects) {
    const s = r.start || r.x || 1, e = r.end || r.to || s;
    const val = valueFn(r);
    for (let col = s; col <= e; col++) {
      const pos = colToRes[col];
      if (pos && (!(pos in map) || val > map[pos])) map[pos] = val;
    }
  }
  return map;
}

function buildSsBfactorMaps() {
  if (!DATA) return null;
  const qaln = DATA.qaln || '', taln = DATA.taln || '';
  const colToA = alnColToResMap(qaln), colToB = alnColToResMap(taln);
  const ssVal = r => r.ss_type === 'helix' ? 2 : r.ss_type === 'strand' ? 1 : 0;
  return {
    A: rectsToResMap(DATA.ssA_alnRects || [], colToA, ssVal),
    B: rectsToResMap(DATA.ssB_alnRects || [], colToB, ssVal),
  };
}

function buildCavityBfactorMaps() {
  if (!DATA) return null;
  const qaln = DATA.qaln || '', taln = DATA.taln || '';
  const colToA = alnColToResMap(qaln), colToB = alnColToResMap(taln);
  const drugVal = r => {
    const d = (r.druggability || '').toLowerCase();
    return d === 'strong' ? 3 : d === 'medium' ? 2 : d === 'weak' ? 1 : 0;
  };
  return {
    A: rectsToResMap(DATA.cavA_alnRects || [], colToA, drugVal),
    B: rectsToResMap(DATA.cavB_alnRects || [], colToB, drugVal),
  };
}

function buildDrugclipBfactorMaps() {
  if (!DATA) return null;
  const qaln = DATA.qaln || '', taln = DATA.taln || '';
  const colToA = alnColToResMap(qaln), colToB = alnColToResMap(taln);
  return {
    A: rectsToResMap(DATA.dcA_alnRects || [], colToA, () => 1),
    B: rectsToResMap(DATA.dcB_alnRects || [], colToB, () => 1),
  };
}

function buildPlmaBfactorMaps() {
  if (!PLMA_DATA || !DATA) return null;
  const mapA = {}, mapB = {};
  const catVal = { specific_a:5, specific_b:5, a_with_family:4, b_with_family:4,
                   pair_exclusive:3, shared_with_family:2, family_only:1 };
  const seqA = PLMA_DATA.gene_a_seq, seqB = PLMA_DATA.gene_b_seq;
  for (const block of (PLMA_DATA.blocks || [])) {
    const val = catVal[block.category] || 0;
    const pA = block.positions?.[seqA];
    if (pA) { for (let r = pA.start; r <= pA.end; r++) { if (!(r in mapA) || val > mapA[r]) mapA[r] = val; } }
    const pB = block.positions?.[seqB];
    if (pB) { for (let r = pB.start; r <= pB.end; r++) { if (!(r in mapB) || val > mapB[r]) mapB[r] = val; } }
  }
  return { A: mapA, B: mapB };
}

function getColoredPdb(mode) {
  if (_pdbBfCache[mode]) return _pdbBfCache[mode];
  const builders = {
    plddt:    [buildPlddtBfactorMaps,    50],
    am:       [buildAmBfactorMaps,       50],
    dam:      [buildDamBfactorMaps,       0],
    aligned:  [buildAlignedBfactorMaps,   0],
    domains:  [buildDomainBfactorMaps,    0],
    ss:       [buildSsBfactorMaps,        0],
    cavities: [buildCavityBfactorMaps,    0],
    drugclip: [buildDrugclipBfactorMaps,  0],
    plma:     [buildPlmaBfactorMaps,      0],
  };
  const entry = builders[mode];
  if (!entry) return null;
  const maps = entry[0]();
  if (!maps) return null;
  const pdb = modifyPdbBfactors(maps.A, maps.B, entry[1]);
  if (pdb) _pdbBfCache[mode] = pdb;
  return pdb;
}

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
    // AlphaMissense: B-factors scaled 0-100 (AM×100)
    // high→low (reverse): red, orange, grey, light grey
    return {
      name: 'uncertainty',
      params: {
        domain: [0, 100],
        list: { kind: 'set', colors: [0xd62728, 0xff7d45, 0xbbbbbb, 0xdddddd] }
      }
    };
  }
  if (m === 'dam') {
    // Delta AM: B-factors |deltaAM|×100, 0-100 range
    return {
      name: 'uncertainty',
      params: {
        domain: [0, 100],
        list: { kind: 'set', colors: [0xd62728, 0xff7d45, 0xbbbbbb, 0xdddddd] }
      }
    };
  }
  if (m === 'aligned') {
    // 4-way: radical=1(green), conservative=0.5(teal), identical=0.25(grey), gap=0(burgundy)
    // Reverse theme (high B → first color): [green, teal, grey, burgundy]
    return {
      name: 'uncertainty',
      params: {
        domain: [0, 1],
        list: { kind: 'set', colors: [0x4CAF50, 0x009688, 0xAAAAAA, 0x8B1A1A] }
      }
    };
  }
  if (m === 'domains') {
    // Distinct color per domain type, 0 = no domain (grey)
    const names = window._domainColorNames || {};
    const nDomains = Object.keys(names).length || 1;
    // Build color list: index 0=grey (no domain), then one color per domain from palette
    // Molstar uncertainty theme with reverse: high→first color
    // Domain [0, nDomains]: 0=grey, 1..N=palette colors
    const colors = [0xcccccc]; // index 0: no domain
    for (let i = 0; i < nDomains; i++) {
      const hex = DOMAIN_PALETTE[i % DOMAIN_PALETTE.length];
      colors.push(parseInt(hex.replace('#', ''), 16));
    }
    // Reverse for uncertainty theme (high B-factor → first color)
    colors.reverse();
    return {
      name: 'uncertainty',
      params: {
        domain: [0, nDomains],
        list: { kind: 'set', colors }
      }
    };
  }
  if (m === 'ss') {
    // Secondary structure: helix=2 (red-pink), strand=1 (yellow), coil=0 (grey)
    return {
      name: 'uncertainty',
      params: {
        domain: [0, 2],
        list: { kind: 'set', colors: [0xFF0066, 0xFFCC00, 0xdddddd] }
      }
    };
  }
  if (m === 'cavities') {
    // Cavity druggability: strong=3, medium=2, weak=1, none=0
    return {
      name: 'uncertainty',
      params: {
        domain: [0, 3],
        list: { kind: 'set', colors: [0xe65100, 0xff9800, 0xffc107, 0xdddddd] }
      }
    };
  }
  if (m === 'drugclip') {
    // DrugCLIP pocket: hit=1 (red), none=0 (white)
    return {
      name: 'uncertainty',
      params: {
        domain: [0, 1],
        list: { kind: 'set', colors: [0xc62828, 0xf7f7f7] }
      }
    };
  }
  if (m === 'plma') {
    // PLMA categories: specific=5, +family=4, pair_excl=3, shared=2, family_only=1, none=0
    // Colors match Nightingale track colors (high→low for reverse theme)
    return {
      name: 'uncertainty',
      params: {
        domain: [0, 5],
        list: { kind: 'set', colors: [0xEF5350, 0xFFA726, 0x26A69A, 0xFFCA28, 0xBDBDBD, 0xEEEEEE] }
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

  try {
    if (mode !== 'uniform') {
      let pdb = getColoredPdb(mode);
      if (pdb) {
        // Filter out hidden chain(s) at PDB level
        pdb = filterPdbChains(pdb, chainVisible.A, chainVisible.B);
        await reloadViewerWith(pdb, true);
        const theme = themeForColorMode(mode);
        await applyColorTheme(theme);
        await renderSelections();
        updateColorLegend(mode);
        return;
      }
      // Fallback to uniform if mode data not available
      currentColorMode = 'uniform';
      const select = document.getElementById('colorBy');
      if (select) select.value = 'uniform';
    }

    if (!plugin || !structureReady) return;
    // For uniform mode, also respect chain visibility
    await applyChainVisibility();
    const theme = themeForColorMode('uniform');
    await applyColorTheme(theme);
    updateColorLegend('uniform');
  } catch(e) {
    console.error('Error applying color theme:', e);
  }
}

// ===== Color Legend =====
function updateColorLegend(mode) {
  const el = document.getElementById('colorLegend');
  if (!el) return;
  if (!mode || mode === 'uniform') { el.style.display = 'none'; return; }

  const legends = {
    plddt: {
      title: 'pLDDT Confidence',
      items: [
        {color:'#0053d6',label:'>90 Very High'},
        {color:'#65cbf3',label:'70-90 Confident'},
        {color:'#ffdb13',label:'50-70 Low'},
        {color:'#ff7d45',label:'≤50 Very Low'},
      ],
      stats: () => {
        const counts = {'> 90':0,'70-90':0,'50-70':0,'≤ 50':0};
        for (const arr of [DATA.plddtA||[], DATA.plddtB||[]]) {
          for (const v of arr) {
            if (v > 90) counts['> 90']++;
            else if (v > 70) counts['70-90']++;
            else if (v > 50) counts['50-70']++;
            else counts['≤ 50']++;
          }
        }
        const total = (DATA.plddtA||[]).length + (DATA.plddtB||[]).length;
        if (!total) return '';
        return Object.entries(counts).map(([k,v])=>`${k}: ${v} (${(v/total*100).toFixed(0)}%)`).join(' · ');
      }
    },
    am: {
      title: 'AlphaMissense Pathogenicity',
      items: [
        {color:'#d62728',label:'Pathogenic (>0.564)'},
        {color:'#ff7d45',label:'Ambiguous (0.34-0.564)'},
        {color:'#bbbbbb',label:'Benign (<0.34)'},
      ],
      stats: () => {
        const counts = {path:0,amb:0,ben:0};
        for (const arr of [DATA.bfactorsA||[], DATA.bfactorsB||[]]) {
          for (const v of arr) {
            if (v > 0.564) counts.path++;
            else if (v > 0.34) counts.amb++;
            else counts.ben++;
          }
        }
        const total = (DATA.bfactorsA||[]).length + (DATA.bfactorsB||[]).length;
        if (!total) return '';
        return `Pathogenic: ${counts.path} (${(counts.path/total*100).toFixed(0)}%) · Ambiguous: ${counts.amb} (${(counts.amb/total*100).toFixed(0)}%) · Benign: ${counts.ben} (${(counts.ben/total*100).toFixed(0)}%)`;
      }
    },
    dam: {
      title: 'Δ AlphaMissense (|A-B| at aligned positions)',
      items: [
        {color:'#d62728',label:'High Δ'},
        {color:'#ff7d45',label:'Medium Δ'},
        {color:'#bbbbbb',label:'Low Δ'},
        {color:'#dddddd',label:'No Δ / Gap'},
      ]
    },
    aligned: {
      title: 'Substitution Type',
      items: [
        {color:'#4CAF50',label:'Radical (diff class)'},
        {color:'#009688',label:'Conservative (same class)'},
        {color:'#AAAAAA',label:'Identical'},
        {color:'#8B1A1A',label:'Gap'},
      ]
    },
    domains: {
      title: 'Domain Regions',
      items: (() => {
        const names = window._domainColorNames || {};
        const items = [];
        for (const [name, idx] of Object.entries(names)) {
          items.push({color: DOMAIN_PALETTE[(idx-1) % DOMAIN_PALETTE.length], label: name});
        }
        items.push({color:'#cccccc', label:'No domain'});
        return items.length > 1 ? items : [{color:'#2ca02c',label:'Domain'},{color:'#cccccc',label:'No domain'}];
      })()
    },
    ss: {
      title: 'Secondary Structure',
      items: [
        {color:'#FF0066',label:'α-Helix'},
        {color:'#FFCC00',label:'β-Strand'},
        {color:'#dddddd',label:'Coil/Loop'},
      ]
    },
    cavities: {
      title: 'Cavity Druggability',
      items: [
        {color:'#e65100',label:'Strong'},
        {color:'#ff9800',label:'Medium'},
        {color:'#ffc107',label:'Weak'},
        {color:'#dddddd',label:'No cavity'},
      ]
    },
    drugclip: {
      title: 'DrugCLIP Binding Pockets',
      items: [
        {color:'#c62828',label:'Pocket hit'},
        {color:'#f7f7f7',label:'No pocket'},
      ]
    },
    plma: {
      title: 'PLMA Conservation Categories',
      items: [
        {color:'#EF5350',label:'Specific'},
        {color:'#FFA726',label:'+ Family'},
        {color:'#26A69A',label:'Pair exclusive'},
        {color:'#FFCA28',label:'Shared w/ family'},
        {color:'#BDBDBD',label:'Family only'},
        {color:'#EEEEEE',label:'None'},
      ]
    },
  };

  const cfg = legends[mode];
  if (!cfg) { el.style.display = 'none'; return; }

  let html = `<strong>${cfg.title}</strong><span style="margin-left:12px">`;
  for (const it of cfg.items) {
    html += `<span style="display:inline-flex;align-items:center;margin-right:10px"><span style="display:inline-block;width:12px;height:12px;border-radius:2px;background:${it.color};margin-right:3px"></span>${it.label}</span>`;
  }
  html += '</span>';
  if (cfg.stats) {
    const s = cfg.stats();
    if (s) html += `<div style="margin-top:4px;color:#666">${s}</div>`;
  }
  el.innerHTML = html;
  el.style.display = '';
}

// ===== Highlight domain dropdowns (checkbox menus) =====
function populateHighlightDropdowns() {
  const wrapA = document.getElementById('hlDomainsA');
  const wrapB = document.getElementById('hlDomainsB');
  if (!wrapA || !wrapB || !DATA) return;

  const lblA = document.getElementById('hlLabelA');
  const lblB = document.getElementById('hlLabelB');
  if (lblA) lblA.textContent = DATA.g1 || 'A';
  if (lblB) lblB.textContent = DATA.g2 || 'B';

  function buildMenu(wrap, domains, chain) {
    wrap.innerHTML = '';
    if (!domains || !domains.length) { wrap.style.display = 'none'; return; }
    const btn = document.createElement('button');
    btn.className = 'hl-menu-btn';
    btn.type = 'button';
    btn.textContent = 'Select domains...';
    btn.style.cssText = 'font-size:inherit;padding:4px 10px;cursor:pointer;min-width:170px;text-align:left;background:#fff;border:1px solid #ccc;border-radius:4px';
    const menu = document.createElement('div');
    menu.className = 'hl-menu-list';
    menu.style.cssText = 'display:none;position:absolute;z-index:100;background:#fff;border:1px solid #ccc;border-radius:4px;box-shadow:0 2px 8px rgba(0,0,0,.15);max-height:280px;overflow-y:auto;min-width:220px;padding:4px 0';
    const seen = new Set();
    for (const d of domains) {
      if (!d.uid) continue;
      const key = `${d.label||d.name}_${d.start}-${d.end}`;
      if (seen.has(key)) continue;
      seen.add(key);
      const item = document.createElement('label');
      item.style.cssText = 'display:flex;align-items:center;gap:6px;padding:5px 12px;cursor:pointer;font-size:13px;white-space:nowrap';
      item.addEventListener('mouseenter', () => { item.style.background='#f0f0f0'; });
      item.addEventListener('mouseleave', () => { item.style.background=''; });
      const cb = document.createElement('input');
      cb.type = 'checkbox';
      cb.value = d.uid;
      cb.dataset.chain = chain;
      cb.style.cssText = 'margin:0;cursor:pointer';
      const swatch = document.createElement('span');
      swatch.style.cssText = `display:inline-block;width:10px;height:10px;border-radius:2px;background:${d.color||'#ccc'};flex-shrink:0`;
      const txt = document.createElement('span');
      txt.textContent = `${d.label||d.name||d.type} (${d.start}-${d.end})`;
      item.appendChild(cb);
      item.appendChild(swatch);
      item.appendChild(txt);
      cb.addEventListener('change', () => {
        const domMap = {};
        for (const dd of domains) { if (dd.uid) domMap[dd.uid] = dd; }
        const dom = domMap[cb.value];
        if (dom) toggleFeature(dom, chain);
        updateMenuBtnLabel(btn, menu);
      });
      menu.appendChild(item);
    }
    wrap.style.position = 'relative';
    wrap.appendChild(btn);
    wrap.appendChild(menu);
    btn.addEventListener('click', (e) => {
      e.stopPropagation();
      const open = menu.style.display !== 'none';
      menu.style.display = open ? 'none' : 'block';
    });
    document.addEventListener('click', (e) => {
      if (!wrap.contains(e.target)) menu.style.display = 'none';
    });
  }
  buildMenu(wrapA, DATA.domainsA, chainIdA);
  buildMenu(wrapB, DATA.domainsB, chainIdB);
}

function updateMenuBtnLabel(btn, menu) {
  const checked = menu.querySelectorAll('input:checked');
  if (checked.length === 0) btn.textContent = 'Select domains...';
  else btn.textContent = `${checked.length} domain${checked.length>1?'s':''} selected`;
}

function syncHighlightDropdowns() {
  for (const [wrapId, chain] of [['hlDomainsA', chainIdA], ['hlDomainsB', chainIdB]]) {
    const wrap = document.getElementById(wrapId);
    if (!wrap) continue;
    const cbs = wrap.querySelectorAll('input[type="checkbox"]');
    for (const cb of cbs) {
      cb.checked = selection.has(selectionKey(chain, cb.value));
    }
    const btn = wrap.querySelector('.hl-menu-btn');
    const menu = wrap.querySelector('.hl-menu-list');
    if (btn && menu) updateMenuBtnLabel(btn, menu);
  }
}

// ===== Alignment selector dropdown =====
function populateAlnSelector() {
  const sel = document.getElementById('alnSelector');
  if (!sel || !DATA) return;
  sel.innerHTML = '<option value="full">Full alignment</option>';
  if (Array.isArray(DATA.domPairs)) {
    DATA.domPairs.forEach((r, i) => {
      const opt = document.createElement('option');
      opt.value = String(i);
      opt.textContent = `${r.Aname} ${r.Arng} × ${r.Bname} ${r.Brng}`;
      sel.appendChild(opt);
    });
  }
  sel.addEventListener('change', async () => {
    if (sel.value === 'full') {
      document.getElementById('backFull').click();
    } else {
      const idx = parseInt(sel.value, 10);
      const r = DATA.domPairs[idx];
      if (!r) return;
      await reloadViewerWith(r.pdb64, true);
      // Reapply color theme (domain PDB may have different B-factors, but theme is preserved)
      if (currentColorMode && currentColorMode !== 'uniform') {
        const theme = themeForColorMode(currentColorMode);
        await applyColorTheme(theme);
      }
      setTmScoreDisplay(r.tm);
      document.getElementById('contextTitle').textContent = `Domain: ${r.Aname} ${r.Arng} × ${r.Bname} ${r.Brng}`;
      selection.clear();
      const domA = (DATA.domainsA||[]).find(d => (d.label===r.Aname) || (`${d.start}-${d.end}`===r.Arng));
      const domB = (DATA.domainsB||[]).find(d => (d.label===r.Bname) || (`${d.start}-${d.end}`===r.Brng));
      if (domA) selection.set(selectionKey(chainIdA, domA.uid), {
        id: domA.uid, chain: chainIdA,
        start: parseInt(domA.start,10), end: parseInt(domA.end,10),
        color: domA.color||'#ffdb13', name: domA.label||domA.name
      });
      if (domB) selection.set(selectionKey(chainIdB, domB.uid), {
        id: domB.uid, chain: chainIdB,
        start: parseInt(domB.start,10), end: parseInt(domB.end,10),
        color: domB.color||'#ff7d45', name: domB.label||domB.name
      });
      await renderSelections();
      syncHighlightDropdowns();
    }
  });
}

// ===== Manual residue selection =====
function setupManualSelection() {
  const inputA = document.getElementById('manualSelA');
  const inputB = document.getElementById('manualSelB');
  const clearBtn = document.getElementById('clearAllHighlights');
  const lblA = document.getElementById('selLabelA');
  const lblB = document.getElementById('selLabelB');
  if (lblA) lblA.textContent = DATA?.g1 || 'A';
  if (lblB) lblB.textContent = DATA?.g2 || 'B';

  function parseRanges(text) {
    const residues = new Set();
    for (const part of text.split(',')) {
      const t = part.trim();
      if (!t) continue;
      const m = t.match(/^(\d+)\s*-\s*(\d+)$/);
      if (m) { for (let i = +m[1]; i <= +m[2]; i++) residues.add(i); }
      else if (/^\d+$/.test(t)) residues.add(+t);
    }
    return residues;
  }

  function applyManual(input, chain) {
    const residues = parseRanges(input.value);
    // Remove existing manual selections for this chain
    for (const [key] of selection) {
      if (key.startsWith(`manual_${chain}_`)) selection.delete(key);
    }
    if (residues.size === 0) { renderSelections(); syncHighlightDropdowns(); return; }
    // Group into contiguous ranges
    const sorted = [...residues].sort((a,b)=>a-b);
    let start = sorted[0], prev = sorted[0];
    for (let i = 1; i <= sorted.length; i++) {
      const curr = i < sorted.length ? sorted[i] : null;
      if (curr === null || curr !== prev + 1) {
        const key = `manual_${chain}_${start}-${prev}`;
        selection.set(key, {
          id: key, chain, start, end: prev,
          color: '#9c27b0', name: `${start}-${prev}`
        });
        if (curr !== null) start = curr;
      }
      if (curr !== null) prev = curr;
    }
    renderSelections();
    syncHighlightDropdowns();
  }

  if (inputA) {
    inputA.addEventListener('keydown', (e) => { if (e.key === 'Enter') applyManual(inputA, chainIdA); });
    inputA.addEventListener('blur', () => applyManual(inputA, chainIdA));
  }
  if (inputB) {
    inputB.addEventListener('keydown', (e) => { if (e.key === 'Enter') applyManual(inputB, chainIdB); });
    inputB.addEventListener('blur', () => applyManual(inputB, chainIdB));
  }
  if (clearBtn) {
    clearBtn.addEventListener('click', () => {
      selection.clear();
      if (inputA) inputA.value = '';
      if (inputB) inputB.value = '';
      renderSelections();
      syncHighlightDropdowns();
    });
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

  // Defer canvas-heavy sections until they're near the viewport
  registerLazySection('similaritySearchSection', initSimilaritySearchSection);
  registerLazySection('familyFeaturesSection', initFamilyFeaturesSection);

  (DATA.domainsA||[]).forEach(d => { if (d.uid) domByUidA[d.uid] = d; });
  (DATA.domainsB||[]).forEach(d => { if (d.uid) domByUidB[d.uid] = d; });
  // Assign distinct colors to domains (same name = same color across A and B)
  assignDomainColors(DATA.domainsA, DATA.domainsB);
  // Re-color server-generated rects with the new palette
  function recolorRects(rects, domains) {
    if (!rects || !domains) return;
    const uidToColor = {};
    for (const d of domains) { if (d.uid && d.color) uidToColor[d.uid] = d.color; }
    for (const r of rects) { if (r.id && uidToColor[r.id]) r.color = uidToColor[r.id]; }
  }
  recolorRects(DATA.domA_alnRects, DATA.domainsA);
  recolorRects(DATA.domB_alnRects, DATA.domainsB);
  recolorRects(DATA.tedA_alnRects, DATA.domainsA);
  recolorRects(DATA.tedB_alnRects, DATA.domainsB);

  setTmScoreDisplay(DATA.tm);
  document.getElementById('contextTitle').textContent = `Full: ${DATA.g1} × ${DATA.g2}`;
  await reloadViewerWith(PDB64_FULL);

  // Defer Nightingale track initialization until alignment section is near-viewport
  registerLazySection('alignmentSection', () => {
    buildSeq();
    setupTrackGroupToggles(); // re-wire chips with actual trackGroupRows
  });
  fillDomainTables();
  await fillDomPairs();
  fillDrugHits();
  populateHighlightDropdowns();
  populateAlnSelector();
  setupManualSelection();

  setupPdbeCollapse();
  setupPdbeControls();
  setupAllCollapsibleSections();
  applyDefaultCollapseStates();
  updateSectionVisibility();

  // Band toggle headers in regions section
  document.querySelectorAll('.band-header[data-band]').forEach(header => {
    header.addEventListener('click', () => {
      const bandId = header.dataset.band;
      const body = document.querySelector(`[data-band-body="${bandId}"]`);
      if (body) {
        const collapsed = body.classList.toggle('collapsed');
        header.classList.toggle('collapsed', collapsed);
      }
    }, { passive: true });
  });

  document.getElementById('colorBy').addEventListener('change', (e)=>colorBy(e.target.value), {passive:true});
  document.getElementById('center').addEventListener('click', ()=>{ if(structureReady){ plugin.canvas3d?.requestCameraReset(); }}, {passive:true});
  document.getElementById('lockViewer').addEventListener('click', ()=>{ toggleViewerLock(); }, {passive:true});
  document.getElementById('focusSelection').addEventListener('click', ()=>{ focusMainViewerSelection(); }, {passive:true});

  // Set initial lock state (delayed to ensure DOM is ready)
  setTimeout(() => {
    setViewerLocked(true);
  }, 100);
  document.getElementById('backFull').addEventListener('click', async ()=>{
    selection.clear();
    pendingHighlightLoci = null;
    pendingSelectionLoci = null;
    // applyChainVisibility re-applies color theme correctly
    await applyChainVisibility();
    setTmScoreDisplay(DATA.tm);
    document.getElementById('contextTitle').textContent = `Full: ${DATA.g1} × ${DATA.g2}`;
    Object.values(trackRefs).forEach(track => {
      if (track && track._originalData) track.data = [...track._originalData];
    });
    await initializeHighlightColors();
    setupHoverInterception();
    await renderSelections();
    syncHighlightDropdowns();
    const alnSel = document.getElementById('alnSelector');
    if (alnSel) alnSel.value = 'full';
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
  
  populateKeyFindingsBanner();
  setupStickyObserver();
  setupStickyMinimize();
  console.log('Report viewer initialized');
}

// Populate the key findings banner with headline stats
function populateKeyFindingsBanner() {
  const banner = document.getElementById('keyFindingsBanner');
  if (!banner) return;

  const conservation = SUMMARY.conservation || {};
  const gene1 = SUMMARY.gene1 || {};
  const gene2 = SUMMARY.gene2 || {};
  const pair = SUMMARY.pair || {};
  const slFunc = SUMMARY.sl_functional || {};

  // Seq identity
  const seqIdent = conservation.min_sequence_identity;
  const kfSeqIdent = document.getElementById('kfSeqIdent');
  if (kfSeqIdent && seqIdent && typeof seqIdent.value === 'number') {
    kfSeqIdent.textContent = (seqIdent.value * 100).toFixed(0) + '%';
    kfSeqIdent.className = 'kf-value ' + (seqIdent.value >= 0.5 ? 'positive' : seqIdent.value >= 0.3 ? 'neutral' : 'negative');
  } else if (kfSeqIdent && DATA.fident != null) {
    kfSeqIdent.textContent = DATA.fident.toFixed(0) + '%';
    kfSeqIdent.className = 'kf-value ' + (DATA.fident >= 50 ? 'positive' : DATA.fident >= 30 ? 'neutral' : 'negative');
  }

  // TM-score
  const kfTm = document.getElementById('kfTmScore');
  if (kfTm && DATA.tm != null) {
    kfTm.textContent = DATA.tm.toFixed(2);
    kfTm.className = 'kf-value ' + (DATA.tm >= 0.5 ? 'positive' : DATA.tm >= 0.3 ? 'neutral' : 'negative');
  }

  // Family size
  const kfFamily = document.getElementById('kfFamilySize');
  if (kfFamily && pair.family_size != null) {
    kfFamily.textContent = pair.family_size + ' members';
  }

  // SL status
  const kfSl = document.getElementById('kfSlStatus');
  if (kfSl) {
    if (slFunc.is_sl === true) { kfSl.textContent = 'Yes'; kfSl.className = 'kf-value positive'; }
    else if (slFunc.is_sl === false) { kfSl.textContent = 'No'; kfSl.className = 'kf-value negative'; }
    else { kfSl.textContent = '?'; kfSl.className = 'kf-value'; }
  }

  // Pockets (cavities + DrugCLIP)
  const kfDrug = document.getElementById('kfDrug');
  if (kfDrug) {
    const nCav = (DATA.cavA_alnRects || []).length + (DATA.cavB_alnRects || []).length;
    const nDc = (DATA.domainsA || []).filter(d => d.type === 'DrugCLIP' || d.raw_type === 'DrugCLIP').length +
                (DATA.domainsB || []).filter(d => d.type === 'DrugCLIP' || d.raw_type === 'DrugCLIP').length;
    const total = nCav + nDc;
    kfDrug.textContent = total;
    kfDrug.className = 'kf-value ' + (total > 0 ? 'positive' : '');
  }

  // Known drugs
  const kfKnownDrugs = document.getElementById('kfKnownDrugs');
  if (kfKnownDrugs) {
    const nDrugs = (gene1.known_drugs || []).length + (gene2.known_drugs || []).length;
    kfKnownDrugs.textContent = nDrugs;
    kfKnownDrugs.className = 'kf-value ' + (nDrugs > 0 ? 'positive' : '');
  }

  banner.style.display = '';
}

// Wire track group toggle chips
function setupTrackGroupToggles() {
  const row = document.getElementById('trackToggleRow');
  if (!row) return;
  const chips = row.querySelectorAll('.track-toggle-chip');
  chips.forEach(chip => {
    const group = chip.dataset.group;
    // Set initial active state from trackGroupState
    chip.classList.toggle('active', !!trackGroupState[group]);
    chip.addEventListener('click', () => {
      if (group === 'default') return; // default is always on
      trackGroupState[group] = !trackGroupState[group];
      chip.classList.toggle('active', trackGroupState[group]);
      const rows = trackGroupRows[group] || [];
      rows.forEach(r => { r.style.display = trackGroupState[group] ? '' : 'none'; });
    }, { passive: true });
  });
}

// Sticky observer for structure controls
function setupStickyObserver() {
  const sticky = document.getElementById('structureStickyControls');
  if (!sticky) return;
  // Insert a sentinel div right before the sticky element
  const sentinel = document.createElement('div');
  sentinel.style.height = '1px';
  sentinel.style.marginBottom = '-1px';
  sticky.parentElement.insertBefore(sentinel, sticky);
  const observer = new IntersectionObserver(([e]) => {
    sticky.classList.toggle('is-stuck', !e.isIntersecting);
  }, { threshold: [1], rootMargin: '-43px 0px 0px 0px' });
  observer.observe(sentinel);
}

// Helper: update TM-score in both main and mini views
function setTmScoreDisplay(val) {
  const txt = val != null ? val.toFixed(3) : '–';
  const el = document.getElementById('tmScore');
  const mini = document.getElementById('tmScoreMini');
  if (el) el.textContent = txt;
  if (mini) mini.textContent = txt;
}

// Sticky minimize / expand wiring
function setupStickyMinimize() {
  const sticky = document.getElementById('structureStickyControls');
  const full = document.getElementById('stickyFull');
  const mini = document.getElementById('stickyMini');
  const btnMin = document.getElementById('stickyMinimize');
  const btnExp = document.getElementById('stickyExpand');
  if (!sticky || !full || !mini || !btnMin || !btnExp) return;

  btnMin.addEventListener('click', () => {
    full.style.display = 'none';
    mini.style.display = '';
    sticky.classList.add('is-mini');
  }, { passive: true });

  btnExp.addEventListener('click', () => {
    mini.style.display = 'none';
    full.style.display = '';
    sticky.classList.remove('is-mini');
  }, { passive: true });
}

