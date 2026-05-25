# DESIGN.md
# Research-Grade PCR Primer Design System
## Full Bioinformatics Pipeline — Scientific & Clinical Grade

**Version:** 1.0.0  
**Status:** Production Design Specification  
**Target Environment:** Python 3.12+ | Linux/macOS/Windows  
**Intended Use:** High-risk biological research, diagnostic primer design, genomic analysis pipelines

---

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Core Design Philosophy](#2-core-design-philosophy)
3. [System Architecture](#3-system-architecture)
4. [Module Specifications](#4-module-specifications)
   - 4.1 [Sequence Validator](#41-sequence-validator)
   - 4.2 [BLAST Engine](#42-blast-engine)
   - 4.3 [MSA Engine](#43-msa-engine)
   - 4.4 [Conserved Region Detector](#44-conserved-region-detector)
   - 4.5 [Primer Generator](#45-primer-generator)
   - 4.6 [Thermodynamics Engine](#46-thermodynamics-engine)
   - 4.7 [Secondary Structure Analyzer](#47-secondary-structure-analyzer)
   - 4.8 [Specificity Checker](#48-specificity-checker)
   - 4.9 [Scoring Engine](#49-scoring-engine)
   - 4.10 [Export Manager](#410-export-manager)
5. [Global MSA Constraint](#5-global-msa-constraint-critical)
6. [Data Flow & Traceability](#6-data-flow--traceability)
7. [Thermodynamic Formulas](#7-thermodynamic-formulas)
8. [Primer Design Rules](#8-primer-design-rules)
9. [Error Handling Strategy](#9-error-handling-strategy)
10. [Logging System](#10-logging-system)
11. [Output Specification](#11-output-specification)
12. [Testing Strategy](#12-testing-strategy)
13. [Performance Design](#13-performance-design)
14. [User Interfaces](#14-user-interfaces)
15. [Security & Reliability](#15-security--reliability)
16. [Project File Structure](#16-project-file-structure)
17. [Scientific Integrity Rules](#17-scientific-integrity-rules)
18. [Dependency Matrix](#18-dependency-matrix)

---

## 1. Project Overview

This system is a **research-grade, fully transparent bioinformatics pipeline** for automated PCR primer design from user-supplied nucleotide sequences.

### 1.1 Intended Use Cases

| Domain | Application |
|--------|-------------|
| Molecular Biology | Target-specific primer design from reference genomes |
| Clinical Diagnostics | Pathogen detection assay development |
| Pharmaceutical R&D | Gene expression and knockdown studies |
| Epidemiology | Multi-strain conserved region identification |
| Forensic Genomics | Species-specific amplification |
| Veterinary Science | Zoonotic pathogen detection |

### 1.2 Non-Goals

This system does **not**:
- Design probes for hybridization arrays
- Perform de novo sequence assembly
- Substitute for wet-lab empirical validation
- Provide clinical diagnostic clearance (use as a research aid only)

### 1.3 Key Guarantees

- ✅ Every computation step is logged and traceable
- ✅ No silent filtering, rejection, or modification of data
- ✅ Global Multiple Sequence Alignment (MSA) only — no local alignment shortcuts
- ✅ All BLAST hits shown, including rejected ones with reasons
- ✅ All primer candidates shown, including rejected ones with reasons
- ✅ Deterministic: same input always produces same output
- ✅ Reproducibility log generated on every run

---

## 2. Core Design Philosophy

### 2.1 Transparency First

Every transformation in the pipeline must produce visible, human-readable intermediate output:

```
Input Sequence
    ↓ [Validation Report]
Cleaned Sequence
    ↓ [BLAST Full Hit Table]
Homologous Sequences
    ↓ [Full MSA Matrix]
Aligned Sequences
    ↓ [Conservation Scores per Position]
Conserved Regions
    ↓ [All Candidate Primers + Rejection Reasons]
Final Primer Set
    ↓ [Scoring Breakdown]
Ranked Primer Report
```

### 2.2 No Hidden Computation

The following behaviors are **strictly prohibited**:

| Prohibited Behavior | Consequence |
|---------------------|-------------|
| Silent BLAST hit removal | System halt with error log |
| Local alignment substituted for global | Hard exception raised |
| Primer rejected without reason | Validation failure |
| Heuristic used without disclosure | Must appear in log as WARNING |
| Missing intermediate output | Pipeline stops, incomplete run flagged |

### 2.3 Biological Correctness

All computations must respect:
- DNA directionality: 5′ → 3′ always
- Watson-Crick complementarity (A↔T, G↔C)
- Reverse complement for reverse primer generation
- Thermodynamic constraints based on nearest-neighbor model or simplified Tm formulas
- Conservation scoring derived only from MSA column statistics

### 2.4 Determinism

- Random seeds are fixed and logged
- All sorting operations are stable
- BLAST results are cached per query hash
- MSA output is identical across runs with identical input

---

## 3. System Architecture

### 3.1 Pipeline Overview

```
┌──────────────────────────────────────────────────────────┐
│                        USER INPUT                         │
│         (FASTA file / raw string / multi-FASTA)          │
└─────────────────────────┬────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│                  SEQUENCE VALIDATOR                       │
│  • FASTA format check                                     │
│  • IUPAC nucleotide validation                            │
│  • RNA → DNA conversion (if needed)                       │
│  • Composition statistics                                 │
└─────────────────────────┬────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│                    BLASTn ENGINE                          │
│  • Query NCBI nt database                                 │
│  • Full hit table output (ALL hits)                       │
│  • Filtered set with rejection reasons                    │
│  • Download homologous sequences                          │
└─────────────────────────┬────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│               MSA ENGINE (GLOBAL ONLY)                    │
│  • ClustalW or Clustal Omega                              │
│  • Full alignment matrix output                           │
│  • Position-wise conservation table                       │
│  • Consensus sequence generation                          │
└─────────────────────────┬────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│             CONSERVED REGION DETECTOR                     │
│  • Entropy scoring per column                             │
│  • Gap density filtering                                  │
│  • Region boundary detection                              │
│  • Conservation evidence table                            │
└─────────────────────────┬────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│                PRIMER DESIGN ENGINE                       │
│  • Forward primer from 5′ conserved region               │
│  • Reverse primer from 3′ conserved region               │
│  • Step-by-step reverse complement shown                  │
│  • All candidates generated and logged                    │
└─────────────────────────┬────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│             THERMODYNAMICS ENGINE                         │
│  • Tm (basic formula)                                     │
│  • Tm (advanced formula)                                  │
│  • GC content validation                                  │
│  • Annealing temperature estimate                         │
└─────────────────────────┬────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│          SECONDARY STRUCTURE ANALYZER                     │
│  • Hairpin detection                                      │
│  • Self-dimer detection                                   │
│  • Cross-dimer detection                                  │
│  • 3′ complementarity check                               │
│  • ΔG stability estimate                                  │
└─────────────────────────┬────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│              SPECIFICITY CHECKER                          │
│  • BLASTn primer sequences against NCBI                   │
│  • Off-target hit mapping                                 │
│  • Pseudogene binding detection                           │
│  • Multi-locus risk scoring                               │
└─────────────────────────┬────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│              SCORING & RANKING ENGINE                     │
│  • Multi-dimensional score (0–100)                        │
│  • Weighted composite with breakdown                      │
│  • Confidence rating                                      │
│  • Risk classification                                    │
└─────────────────────────┬────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│                 OUTPUT & EXPORT                           │
│  • TXT / CSV / JSON / FASTA / PDF                         │
│  • Full traceability report                               │
│  • Reproducibility log                                    │
└──────────────────────────────────────────────────────────┘
```

### 3.2 Component Interaction

```
pipeline_runner.py
    ├── sequence_validator.py       → ValidatedSequence
    ├── blast_engine.py             → BlastResultSet
    ├── msa_engine.py               → AlignmentMatrix
    ├── conserved_region_detector.py → ConservedRegionList
    ├── primer_generator.py         → PrimerCandidateList
    ├── thermodynamics.py           → ThermoProfile
    ├── secondary_structure.py      → StructureRiskReport
    ├── specificity_checker.py      → SpecificityReport
    ├── scoring_engine.py           → ScoredPrimerSet
    └── export_manager.py           → FinalReportBundle
```

---

## 4. Module Specifications

### 4.1 Sequence Validator

**File:** `src/modules/sequence_validator.py`

#### Responsibilities
- Parse FASTA format (single and multi-FASTA)
- Accept raw nucleotide strings
- Validate against IUPAC nucleotide codes
- Detect and optionally convert RNA sequences (U → T)
- Identify and report ambiguous bases (R, Y, S, W, K, M, B, D, H, V, N)
- Report exact character position of any invalid base

#### Input
```python
@dataclass
class ValidationInput:
    raw_input: str              # FASTA string or raw sequence
    input_format: str           # "fasta" | "raw" | "multi_fasta"
    allow_iupac: bool = True    # Whether to accept ambiguous codes
    convert_rna: bool = True    # Auto-convert U → T
```

#### Output
```python
@dataclass
class ValidationResult:
    is_valid: bool
    cleaned_sequence: str
    sequence_id: str
    sequence_length: int
    gc_content: float               # As percentage (0–100)
    composition: dict[str, int]     # {'A': n, 'T': n, 'G': n, 'C': n, ...}
    ambiguous_bases: list[tuple]    # [(position, char), ...]
    errors: list[str]               # Exact error messages with positions
    warnings: list[str]
    was_rna_converted: bool
    validation_report: str          # Human-readable full report
```

#### Validation Rules

| Rule | Action on Failure |
|------|-------------------|
| Characters not in IUPAC set | Reject: log position and character |
| Sequence length < 50 bp | Warning: may be too short for BLASTn |
| Sequence length < 18 bp | Hard reject: too short for any primer |
| GC% < 20% or > 80% | Warning: extreme composition |
| Contains only one nucleotide type | Warning: low complexity |
| FASTA header missing | Warning: sequence labeled as "unknown" |

---

### 4.2 BLAST Engine

**File:** `src/modules/blast_engine.py`

#### Responsibilities
- Submit BLASTn query to NCBI via Biopython `NCBIWWW.qblast()` or REST API
- Retrieve and parse all hits
- Apply identity/coverage filters with full logging
- Download FASTA sequences for filtered hits
- Cache results to disk by query hash

#### Input
```python
@dataclass
class BlastConfig:
    database: str = "nt"
    program: str = "blastn"
    identity_threshold: float = 95.0
    coverage_threshold: float = 80.0
    max_hits: int = 100
    e_value_cutoff: float = 1e-10
    entrez_query: str = ""          # Optional: filter by organism
    use_cache: bool = True
```

#### Output — Per Hit (ALL hits, including rejected)

```python
@dataclass
class BlastHit:
    accession: str
    organism: str
    alignment_length: int
    query_coverage: float       # Percentage
    identity_percent: float
    mismatches: int
    gaps: int
    e_value: float
    bit_score: float
    aligned_segment_preview: str  # First 80 chars of alignment
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    strand: str                 # "plus/plus" | "plus/minus"
    is_accepted: bool
    rejection_reason: str | None  # MUST be populated if rejected
    raw_score: float
```

#### Filter Logging Example

```
[BLAST FILTER] Hit: NM_001234.5 | Homo sapiens
  Identity: 91.2% < threshold 95.0% → REJECTED
  Coverage: 98.1% ✓
  E-value: 2.3e-45 ✓
  Reason: Identity below threshold (91.2% < 95.0%)
```

#### Critical Rules

- `is_accepted = False` entries must **never** be silently dropped
- Every hit must appear in the final report under "All BLAST Results"
- Duplicate removal must be logged: `"Duplicate: NM_001234.5 already present → skipped"`

---

### 4.3 MSA Engine

**File:** `src/modules/msa_engine.py`

#### ⚠️ CRITICAL CONSTRAINT: GLOBAL ALIGNMENT ONLY

This module must enforce:
- ClustalW or Clustal Omega in **global alignment mode**
- Every sequence aligned from start to end
- No sequence trimming before alignment
- No local alignment fallback under any condition

If a global alignment tool is unavailable, the system must:
1. Stop execution
2. Log a `CRITICAL` error
3. Refuse to proceed with local alignment

#### Input
```python
@dataclass
class MSAInput:
    sequences: list[SeqRecord]      # All BLAST-retrieved sequences + query
    aligner: str = "clustalo"       # "clustalw" | "clustalo"
    output_format: str = "clustal"
```

#### Output — Full Alignment Matrix

```python
@dataclass
class AlignmentMatrix:
    aligned_sequences: list[AlignedSequence]   # All seqs, same length
    alignment_length: int
    num_sequences: int
    conservation_table: list[PositionStats]
    strict_consensus: str
    majority_consensus: str
    weighted_consensus: str
    gap_distribution: list[float]   # Gap % per position
    full_alignment_text: str        # Complete Clustal format text
```

#### Required Console/Log Output Format

**Full Alignment Matrix (every sequence shown):**
```
=== FULL MSA ALIGNMENT (GLOBAL) ===
Total sequences: 12 | Alignment length: 1,847 bp

Seq_Query      : ATGCTAGCTAGCT---AAGTCGATCG
Seq_NM_001234  : ATG-TAGCT-GCTAACAAGTCGATCG
Seq_NM_005678  : ATGCTAG-TAGCTAACAAGTCGATCG
...
(all 12 sequences shown)
```

**Position-wise Conservation Table:**
```
=== CONSERVATION TABLE (ALL POSITIONS) ===
Pos   A%    T%    G%    C%    Gap%  Entropy  Status
001   100   0     0     0     0     0.000    CONSERVED
002   0     0     0     100   0     0.000    CONSERVED
003   0     0     85    15    0     0.610    HIGH_CONS
004   40    60    0     0     5     0.971    MODERATE
...
(all 1847 positions shown)
```

#### Conservation Status Thresholds

| Status | Conservation % | Gap % |
|--------|---------------|-------|
| CONSERVED | 100% | 0% |
| HIGH_CONSERVATION | ≥ 90% | ≤ 5% |
| MODERATE | 70–90% | ≤ 10% |
| VARIABLE | < 70% | any |
| GAP_RICH | any | > 20% |

---

### 4.4 Conserved Region Detector

**File:** `src/modules/conserved_region_detector.py`

#### Selection Criteria

A region qualifies as conserved **only if ALL** of the following are satisfied:

| Criterion | Threshold |
|-----------|-----------|
| Conservation score | ≥ 90% across all sequences |
| Gap frequency | ≤ 5% |
| Minimum length | ≥ 18 bp (for primer extraction) |
| Stability | Consistent across all aligned sequences |

No region may be selected without showing raw evidence from the MSA.

#### Output Per Region

```python
@dataclass
class ConservedRegion:
    region_id: str
    start_position: int         # Alignment index (1-based)
    end_position: int
    length: int
    conservation_score: float
    gap_frequency: float
    entropy_mean: float
    consensus_sequence: str
    per_organism_sequences: dict[str, str]  # {organism: sequence}
    is_suitable_forward: bool   # Suitable for forward primer
    is_suitable_reverse: bool   # Suitable for reverse primer
    evidence_table: str         # Full MSA columns for this region
```

#### Output Log Example

```
=== CONSERVED REGION CR-001 ===
Position: 142–178 (alignment index)
Length: 37 bp
Conservation: 96.4%
Gap frequency: 1.2%
Entropy (mean): 0.083

Sequence per organism:
  Homo sapiens (NM_001234): ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
  Mus musculus (NM_005678): ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
  Rattus norvegicus (NC_00): ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
  ...

Consensus: ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
Status: ACCEPTED → Suitable for FORWARD primer
```

---

### 4.5 Primer Generator

**File:** `src/modules/primer_generator.py`

#### Primer Generation Strategy

- **Forward primer:** Extracted from 5′ conserved region
- **Reverse primer:** Extracted from 3′ conserved region, then reverse-complemented

#### Reverse Complement — Step-by-Step Display (MANDATORY)

```
=== REVERSE PRIMER DERIVATION ===
1. Source region sequence (5′→3′):
   ATGCTAGCTAGCTAGCTAGCT

2. Complement (A↔T, G↔C):
   TACGATCGATCGATCGATCGA

3. Reverse (read 3′→5′ as 5′→3′):
   AGCTAGCTAGCTAGCTAGCAT

4. Final reverse primer (5′→3′):
   AGCTAGCTAGCTAGCTAGCAT
```

#### Per-Candidate Output (ALL candidates, including rejected)

```python
@dataclass
class PrimerCandidate:
    candidate_id: str
    direction: str              # "FORWARD" | "REVERSE"
    sequence: str               # 5′→3′
    source_region_id: str
    source_start: int           # Alignment position
    source_end: int
    length: int
    gc_content: float
    tm_basic: float
    tm_advanced: float
    mismatch_positions: list[int]
    conservation_score: float
    secondary_structure_risk: str   # "LOW" | "MEDIUM" | "HIGH"
    repeat_detected: bool
    repeat_sequence: str | None
    gc_clamp: bool
    is_accepted: bool
    rejection_reasons: list[str]    # MUST be populated if rejected
    organism_contributions: dict[str, str]
```

---

### 4.6 Thermodynamics Engine

**File:** `src/modules/thermodynamics.py`

#### Formula 1 — Basic Wallace Rule

```
Tm = 2(A + T) + 4(G + C)
```

Where A, T, G, C are the counts of each nucleotide in the primer.

#### Formula 2 — Adjusted Formula

```
Tm = 64.9 + (41 × (G + C − 16.4)) / (A + T + G + C)
```

#### Annealing Temperature

```
Ta = Tm − 5°C   (general guideline, disclosed as heuristic)
```

#### GC Content

```
GC% = (G + C) / (A + T + G + C) × 100
```

#### Per-Primer Output

```
=== THERMODYNAMIC ANALYSIS: Primer FWD-001 ===
Sequence: ATGCTAGCTAGCTAGCTAGC
Length: 20 bp

Composition: A=5  T=4  G=6  C=5

Tm (Basic Wallace):    Tm = 2(9) + 4(11) = 18 + 44 = 62°C
Tm (Advanced formula): Tm = 64.9 + (41 × (11 − 16.4)) / 20
                            = 64.9 + (41 × −5.4) / 20
                            = 64.9 − 11.07 = 53.83°C

GC Content: 55.0% → PASS (40–60% range)
Tm difference with pair: 1.2°C → PASS (≤ 2°C)
Annealing Temperature (estimate): 48.83°C
```

---

### 4.7 Secondary Structure Analyzer

**File:** `src/modules/secondary_structure.py`

#### Detection Targets

| Structure | Detection Method |
|-----------|-----------------|
| Hairpin | Self-complementarity sliding window, ΔG estimate |
| Self-dimer | Full self vs self alignment |
| Cross-dimer | Forward vs Reverse alignment |
| 3′ Binding | Last 5 bp complementarity check |

#### ΔG Estimation

Simplified ΔG estimate using nearest-neighbor parameters (disclosed as approximation):

```
ΔG_hairpin ≈ ΔH − T × ΔS   (nearest-neighbor stacking values)
```

Note: If Vienna RNA package or equivalent is unavailable, simplified estimation is used and **must be disclosed** in the report.

#### Output Per Primer

```python
@dataclass
class StructureRiskReport:
    primer_id: str
    hairpin_detected: bool
    hairpin_sequence: str | None
    hairpin_binding_positions: tuple | None
    hairpin_delta_g: float | None
    self_dimer_detected: bool
    self_dimer_binding_positions: list[tuple] | None
    cross_dimer_detected: bool          # vs paired primer
    cross_dimer_binding_positions: list[tuple] | None
    three_prime_binding: bool
    three_prime_complement_seq: str | None
    overall_risk: str                   # "LOW" | "MEDIUM" | "HIGH"
    pass_fail: bool
    failure_explanation: str | None
```

#### Severity Thresholds

| Risk Level | Criteria |
|------------|----------|
| LOW | No hairpin; dimers with ΔG > −1.0 kcal/mol |
| MEDIUM | Hairpin with ΔG > −2.0; weak dimer ΔG > −3.0 |
| HIGH | Hairpin ΔG ≤ −2.0; 3′ self-complementarity ≥ 3 bp |
| FAIL | Any 3′ binding ≥ 4 bp; cross-dimer ΔG ≤ −5.0 |

---

### 4.8 Specificity Checker

**File:** `src/modules/specificity_checker.py`

#### Process

1. Submit each primer sequence to BLASTn (short sequence mode)
2. Identify all genomic binding locations
3. Map off-target hits
4. Flag pseudogene hits
5. Compute risk score

#### Output Per Primer

```python
@dataclass
class SpecificityReport:
    primer_id: str
    total_hits: int
    on_target_hits: list[BlastHit]
    off_target_hits: list[BlastHit]     # ALL shown, not just top N
    pseudogene_hits: list[BlastHit]
    multi_locus_risk: bool
    max_off_target_identity: float
    specificity_score: float            # 0–100
    risk_level: str                     # "LOW" | "MEDIUM" | "HIGH"
    recommendation: str
```

#### Off-Target Risk Classification

| Risk | Condition |
|------|-----------|
| LOW | 1 genomic hit, on-target only |
| MEDIUM | ≤ 3 hits, off-target identity < 80% |
| HIGH | > 3 hits, or any off-target identity ≥ 85% |
| CRITICAL | Pseudogene or repeat region hit with > 90% identity |

---

### 4.9 Scoring Engine

**File:** `src/modules/scoring_engine.py`

#### Scoring Dimensions

| Dimension | Weight | Source |
|-----------|--------|--------|
| Conservation Score | 25% | MSA column statistics |
| Thermodynamic Stability | 20% | Tm, GC%, Ta |
| Specificity Score | 25% | BLASTn primer check |
| Secondary Structure | 15% | Hairpin / dimer ΔG |
| Sequence Quality | 15% | Length, repeats, GC clamp |

#### Composite Score Formula

```
FinalScore = Σ (weight_i × score_i)  for i in dimensions

Where each score_i is normalized to 0–100.
```

#### Per-Primer Score Output

```
=== FINAL SCORE: Primer FWD-001 ===
Conservation Score:       92.0 × 0.25 = 23.00
Thermodynamic Stability:  88.5 × 0.20 = 17.70
Specificity Score:        95.0 × 0.25 = 23.75
Secondary Structure:      85.0 × 0.15 = 12.75
Sequence Quality:         90.0 × 0.15 = 13.50
────────────────────────────────────────────────
COMPOSITE SCORE: 90.70 / 100
Confidence Rating: HIGH
Risk Level: LOW
Recommendation: ACCEPTED — Suitable for experimental use
```

---

### 4.10 Export Manager

**File:** `src/modules/export_manager.py`

#### Output Formats

| Format | Contents | Filename Pattern |
|--------|----------|-----------------|
| TXT | Full human-readable pipeline report | `report_<hash>.txt` |
| CSV | Primer table, BLAST hits, scoring | `primers_<hash>.csv` |
| JSON | Complete structured pipeline data | `pipeline_<hash>.json` |
| FASTA | Final selected primers | `primers_<hash>.fasta` |
| PDF | Scientific report with tables | `report_<hash>.pdf` |
| LOG | Full execution log | `run_<timestamp>.log` |

#### JSON Structure

```json
{
  "run_id": "abc123",
  "timestamp": "2024-01-15T14:32:00Z",
  "input_sequence_hash": "sha256:...",
  "validation": { ... },
  "blast_results": {
    "all_hits": [ ... ],
    "accepted_hits": [ ... ],
    "rejected_hits": [ ... ]
  },
  "msa": {
    "full_alignment": "...",
    "conservation_table": [ ... ],
    "consensus": { ... }
  },
  "conserved_regions": [ ... ],
  "primer_candidates": {
    "all_candidates": [ ... ],
    "accepted": [ ... ],
    "rejected": [ ... ]
  },
  "thermodynamics": [ ... ],
  "secondary_structure": [ ... ],
  "specificity": [ ... ],
  "scores": [ ... ],
  "final_primers": {
    "forward": { ... },
    "reverse": { ... },
    "amplicon_size_bp": 0
  },
  "reproducibility_log": "..."
}
```

---

## 5. Global MSA Constraint (Critical)

### Hard Rules

```
╔══════════════════════════════════════════════════════════════╗
║         MANDATORY: GLOBAL ALIGNMENT ONLY                    ║
╠══════════════════════════════════════════════════════════════╣
║  ✔ All sequences aligned end-to-end                         ║
║  ✔ Every alignment position is accounted for                ║
║  ✔ All gaps explicitly shown as "-"                         ║
║  ✔ No sequence exclusion at any stage                       ║
║  ✗ Local alignment mode is FORBIDDEN                        ║
║  ✗ Heuristic trimming before alignment is FORBIDDEN         ║
║  ✗ Partial alignment output is FORBIDDEN                    ║
╚══════════════════════════════════════════════════════════════╝
```

### Tool Configuration

#### ClustalW Global Mode
```bash
clustalw2 -INFILE=sequences.fasta \
           -TYPE=DNA \
           -ALIGN \
           -OUTPUT=CLUSTAL \
           -OUTFILE=alignment.aln
```

#### Clustal Omega Global Mode
```bash
clustalo --infile sequences.fasta \
         --outfile alignment.aln \
         --outfmt clustal \
         --full \
         --force
```

### Fallback Behavior

If neither ClustalW nor Clustal Omega is found:
1. System raises `AlignmentToolNotFoundError`
2. Logs `CRITICAL: No global alignment tool available`
3. Provides installation instructions
4. **Does not fall back to Biopython PairwiseAligner (local)**
5. Pipeline halts cleanly

---

## 6. Data Flow & Traceability

### 6.1 Traceability Chain

Every data object carries a lineage record:

```python
@dataclass
class DataLineage:
    object_id: str
    source_module: str
    created_at: str         # ISO timestamp
    input_hash: str         # SHA-256 of input data
    parameters_used: dict
    parent_ids: list[str]   # IDs of objects this was derived from
    transformations: list[str]  # List of operations applied
```

### 6.2 Full Traceability Report (Required in Every Run)

```
══════════════════════════════════════════════════════
         FULL TRACEABILITY REPORT
══════════════════════════════════════════════════════
Run ID: run_20240115_143200_abc123
Input hash: sha256:7f3a9b...
Timestamp: 2024-01-15 14:32:00 UTC

[1] VALIDATION
    Input: raw_sequence (1,203 bp)
    Cleaned: 1,203 bp | GC: 52.3%
    Warnings: none

[2] BLAST RESULTS (ALL 47 hits shown)
    Accepted: 12 | Rejected: 35
    Rejection reasons logged: ✓

[3] MSA (Global Clustal Omega)
    Sequences aligned: 13 (12 + query)
    Alignment length: 1,847 bp
    Full matrix: ✓ (13 × 1,847)

[4] CONSERVED REGIONS
    Detected: 8 | Used for primers: 2
    Rejected regions: 6 (reasons logged)

[5] PRIMER CANDIDATES
    Generated: 24 | Accepted: 3 pairs
    Rejected: 21 (reasons logged for all)

[6] THERMODYNAMICS
    All calculations shown step-by-step: ✓

[7] SECONDARY STRUCTURE
    Checked: 6 primers | Passed: 4 | Failed: 2

[8] SPECIFICITY
    BLASTn run for 4 primers
    Off-target hits: all shown

[9] FINAL SELECTION
    Forward primer: FWD-003 (Score: 91.2/100)
    Reverse primer: REV-002 (Score: 89.7/100)
    Amplicon size: 412 bp

REPRODUCIBILITY: All parameters fixed and logged.
══════════════════════════════════════════════════════
```

---

## 7. Thermodynamic Formulas

### 7.1 Basic Wallace Rule (always computed)

```
Tm = 2(A + T) + 4(G + C)
```

- A = count of adenines in primer
- T = count of thymines in primer
- G = count of guanines in primer
- C = count of cytosines in primer

### 7.2 Adjusted Tm Formula (always computed)

```
Tm = 64.9 + [41 × (G + C − 16.4)] / (A + T + G + C)
```

### 7.3 Primer Pair Compatibility

```
|Tm_forward − Tm_reverse| ≤ 2°C   → PASS
|Tm_forward − Tm_reverse| > 2°C   → WARNING (logged)
|Tm_forward − Tm_reverse| > 5°C   → FAIL (primer pair rejected)
```

### 7.4 Disclosure Requirement

Both formulas must be computed for every primer. If results differ significantly:

```
NOTE: Basic Tm = 62°C; Advanced Tm = 53.8°C
Discrepancy of 8.2°C detected.
This is expected for primers near salt/length assumptions.
Advanced formula used for final decision.
```

---

## 8. Primer Design Rules

### 8.1 Mandatory Constraints

| Parameter | Rule | On Violation |
|-----------|------|--------------|
| Length | 18–25 bp | Reject with reason |
| GC content | 40–60% | Reject with actual GC% |
| Tm difference (pair) | ≤ 2°C | Warning at >2°C, reject at >5°C |
| Repeated bases | ≤ 3 consecutive identical | Reject, show repeat |
| Dinucleotide repeats | No ATATAT, CGCGCG, etc. | Reject, show pattern |
| GC clamp at 3′ end | G or C preferred | Warning if absent |
| Self-complementarity | No significant hairpin | Reject if ΔG ≤ −2 kcal/mol |
| Primer-dimer risk | Low cross-dimer ΔG | Reject if ΔG ≤ −5 kcal/mol |

### 8.2 Preferred Properties (Scored, Not Mandatory)

| Property | Preference | Score Impact |
|----------|------------|-------------|
| GC clamp | 1–2 G/C at 3′ end | +5 points |
| Balanced distribution | No 5-base runs | +3 points |
| Tm 55–65°C | Optimal PCR range | +5 points |
| Length 20–22 bp | Optimal length | +3 points |

### 8.3 Rejection Logging (All rejections logged)

```
[PRIMER REJECTED] Candidate FWD-007
  Sequence: ATATATATATATAT GCGCGC
  Reason 1: Dinucleotide repeat detected → ATATAT at positions 1–6
  Reason 2: GC content 64.3% > 60% threshold
  Length: 21 bp (OK)
  Decision: REJECTED (2 violations)
```

---

## 9. Error Handling Strategy

### 9.1 Error Taxonomy

| Error Type | Severity | Action |
|------------|----------|--------|
| Invalid FASTA format | CRITICAL | Stop; show parse error location |
| Invalid nucleotide character | ERROR | Stop; show position |
| BLAST network failure | ERROR | Retry 3× with backoff; use cache if available |
| BLAST returns 0 results | WARNING | Continue with only query sequence |
| Alignment tool not found | CRITICAL | Stop; show install instructions |
| Alignment fails to converge | ERROR | Stop; log alignment parameters |
| No conserved regions found | ERROR | Stop; report coverage statistics |
| No primers pass filters | ERROR | Stop; show all rejections |
| Thermodynamic computation error | ERROR | Stop; show formula inputs |
| Export write failure | WARNING | Try alternate format; notify user |

### 9.2 Exception Classes

```python
class PrimerDesignBaseError(Exception): pass

class SequenceValidationError(PrimerDesignBaseError):
    def __init__(self, message, position=None, character=None): ...

class BLASTFailureError(PrimerDesignBaseError):
    def __init__(self, message, attempt=None, response_code=None): ...

class AlignmentToolNotFoundError(PrimerDesignBaseError):
    def __init__(self, tool_name, install_instructions): ...

class AlignmentFailureError(PrimerDesignBaseError):
    def __init__(self, message, sequences_attempted=None): ...

class NoConservedRegionError(PrimerDesignBaseError):
    def __init__(self, coverage_stats): ...

class NoPrimerCandidatesError(PrimerDesignBaseError):
    def __init__(self, rejection_summary): ...

class ThermodynamicsError(PrimerDesignBaseError):
    def __init__(self, formula, inputs, trace): ...
```

### 9.3 Retry Strategy for API Calls

```python
RETRY_DELAYS = [5, 30, 120]   # seconds (exponential-ish)
MAX_RETRIES = 3

# On all 3 failures: fall back to cached result if available,
# otherwise raise BLASTFailureError with full context.
```

---

## 10. Logging System

### 10.1 Log Levels

| Level | Usage |
|-------|-------|
| DEBUG | Per-position MSA scores, raw alignment output |
| INFO | Module entry/exit, accepted results, statistics |
| WARNING | Threshold borderline cases, approximations used |
| ERROR | Recoverable failures, rejected pipeline stages |
| CRITICAL | Unrecoverable failure, pipeline halt |

### 10.2 Log Entry Format

```
2024-01-15 14:32:00.123 | INFO     | blast_engine:run_blast:89 | BLAST query submitted | query_hash=7f3a9b | db=nt | hits=47
2024-01-15 14:32:45.001 | WARNING  | primer_generator:check_gc:204 | GC content borderline | primer=FWD-003 | gc=59.1%
2024-01-15 14:33:01.872 | ERROR    | msa_engine:align:67 | Clustal Omega not found | PATH=/usr/bin | fallback=NONE
2024-01-15 14:33:01.873 | CRITICAL | pipeline_runner:run:45 | Pipeline halted: global alignment unavailable
```

### 10.3 Mandatory Log Events

- Input sequence hash on every run
- Every BLAST hit (accepted and rejected)
- Full MSA statistics
- Every conserved region (accepted and rejected)
- Every primer candidate (accepted and rejected, with reasons)
- All Tm calculations
- Secondary structure results per primer
- Final selection reasoning

---

## 11. Output Specification

### 11.1 Console Output (always printed)

```
╔══════════════════════════════════════════════════════════════╗
║         PCR PRIMER DESIGN SYSTEM v1.0.0                     ║
║         Research Grade | Full Transparency Mode             ║
╚══════════════════════════════════════════════════════════════╝

[1/9] SEQUENCE VALIDATION ............ ✓ PASS (1,203 bp | GC: 52.3%)
[2/9] BLASTn SEARCH .................. ✓ 47 hits | 12 accepted | 35 rejected
[3/9] MULTIPLE SEQUENCE ALIGNMENT .... ✓ Global | 13 seqs | 1,847 bp
[4/9] CONSERVED REGIONS .............. ✓ 8 found | 2 suitable for primers
[5/9] PRIMER GENERATION .............. ✓ 24 candidates | 6 passed filters
[6/9] THERMODYNAMICS ................. ✓ Tm range: 58.2–62.4°C
[7/9] SECONDARY STRUCTURE ............ ✓ 4/6 primers pass structure check
[8/9] SPECIFICITY CHECK .............. ✓ 4 primers checked | 3 specific
[9/9] FINAL SCORING .................. ✓ Best pair selected

══════════════════════════════════════════════════════════════
SELECTED PRIMERS:
  Forward: ATGCTAGCTAGCTAGCTAGC   Tm: 60.2°C   GC: 55%   Score: 91.2
  Reverse: GCTAGCTAGCTAGCTAGCAT   Tm: 59.6°C   GC: 55%   Score: 89.7
  Amplicon: 412 bp
══════════════════════════════════════════════════════════════
Report saved: report_abc123.pdf
```

### 11.2 Required Report Sections (All Formats)

1. Full MSA alignment (all sequences, all positions)
2. Position-wise conservation table (every position)
3. All conserved regions (accepted and rejected)
4. All primer candidates (accepted and rejected, with reasons)
5. Thermodynamic analysis per primer (step-by-step)
6. Secondary structure analysis per primer
7. Specificity BLASTn results per primer
8. Final scoring breakdown per primer
9. Selected primer pair summary
10. Amplicon size and position
11. Full reproducibility log

---

## 12. Testing Strategy

### 12.1 Test Categories

| Category | Description | Expected Outcome |
|----------|-------------|-----------------|
| valid_standard | Standard 500 bp coding sequence | Full pipeline success |
| valid_fasta | Well-formatted FASTA input | Correct parsing |
| invalid_chars | Sequence with X, Z, 1, spaces | Validation error at exact position |
| rna_input | Sequence with U instead of T | Auto-converted, warning logged |
| too_short | 10 bp sequence | Hard rejection at validation |
| extreme_gc_high | 90% GC content | Warning; pipeline continues |
| extreme_gc_low | 10% GC content | Warning; pipeline continues |
| repetitive_dna | ATATATATATATAT... | Primers rejected; logged |
| blast_failure | Simulated network error | Retry logic triggered; fallback used |
| no_blast_hits | Artificial sequence with no homologs | No-hits handled gracefully |
| alignment_failure | Conflicting sequences | Alignment error logged cleanly |
| no_conserved_region | Fully variable alignment | Error with statistics |
| multi_fasta | 5-sequence FASTA file | All sequences processed |
| large_sequence | 50,000 bp genome segment | Performance within limits |

### 12.2 Test Framework

```bash
pytest tests/ -v --cov=src --cov-report=html
```

- **Coverage target:** ≥ 90% line coverage
- **All tests must be deterministic:** no random variation between runs
- **Integration tests** must use real BLAST/MSA outputs (cached for reproducibility)

---

## 13. Performance Design

### 13.1 Optimization Techniques

| Technique | Application |
|-----------|-------------|
| BLAST result caching | Store results by SHA-256 hash of query sequence |
| MSA result caching | Cache alignment by sequence set hash |
| Multiprocessing | Parallel secondary structure checks per primer |
| Async API calls | Non-blocking BLAST requests |
| Memory mapping | Large sequence files > 10 MB |

### 13.2 Performance Targets

| Operation | Target Time |
|-----------|-------------|
| Sequence validation | < 1 second |
| BLAST search (remote) | < 120 seconds |
| MSA (≤ 50 sequences) | < 30 seconds |
| Primer generation | < 5 seconds |
| Full pipeline (typical) | < 5 minutes |
| Report export (PDF) | < 10 seconds |

---

## 14. User Interfaces

### 14.1 Command-Line Interface

```bash
python -m primer_designer \
    --input sequence.fasta \
    --output-dir ./results \
    --format all \
    --blast-identity 95.0 \
    --blast-coverage 80.0 \
    --min-conservation 0.90 \
    --primer-length-min 18 \
    --primer-length-max 25 \
    --gc-min 40.0 \
    --gc-max 60.0 \
    --verbose
```

**Flags:**
- `--dry-run` — Validate input and show parameters only
- `--no-blast` — Skip BLAST; use provided sequences via `--sequences-file`
- `--offline` — Use cached BLAST results only
- `--report-only` — Reformat a previous run from JSON cache

### 14.2 Streamlit Web App

Pages:
1. **Upload** — FASTA upload, raw sequence input, or example sequence
2. **Validation** — Real-time sequence statistics display
3. **BLAST** — Live BLAST progress, hit table
4. **Alignment** — Interactive MSA viewer (scrollable, color-coded)
5. **Conservation** — Position-by-position conservation chart
6. **Primers** — Candidate table with filter/sort
7. **Thermodynamics** — Tm visualization, pair comparison
8. **Secondary Structure** — Structure diagram per primer
9. **Specificity** — Off-target hit map
10. **Report** — Download links for all formats

---

## 15. Security & Reliability

### 15.1 Input Sanitization

- Maximum input sequence length: 500,000 bp (configurable)
- FASTA header stripped of shell-special characters before any subprocess call
- All external tool calls use argument lists (not shell=True)
- Temporary files written to isolated temp directory, cleaned on exit

### 15.2 API Protection

- NCBI BLAST rate limiting: max 3 requests/second (NCBI policy)
- API key support via environment variable: `NCBI_API_KEY`
- Retry with exponential backoff: 5s, 30s, 120s
- Timeout: 180 seconds per BLAST request

### 15.3 Memory Limits

- Sequences > 500 kb: streamed, not fully loaded into memory
- MSA matrices > 10,000 × 10,000: chunked processing with logging
- Output files: size-checked before write

---

## 16. Project File Structure

```
primer_designer/
├── src/
│   ├── __init__.py
│   ├── pipeline_runner.py          # Orchestrates full pipeline
│   └── modules/
│       ├── __init__.py
│       ├── sequence_validator.py
│       ├── blast_engine.py
│       ├── msa_engine.py
│       ├── conserved_region_detector.py
│       ├── primer_generator.py
│       ├── thermodynamics.py
│       ├── secondary_structure.py
│       ├── specificity_checker.py
│       ├── scoring_engine.py
│       └── export_manager.py
├── tests/
│   ├── __init__.py
│   ├── test_sequence_validator.py
│   ├── test_blast_engine.py
│   ├── test_msa_engine.py
│   ├── test_conserved_region_detector.py
│   ├── test_primer_generator.py
│   ├── test_thermodynamics.py
│   ├── test_secondary_structure.py
│   ├── test_specificity_checker.py
│   ├── test_scoring_engine.py
│   └── fixtures/
│       ├── sample_sequences/
│       ├── mock_blast_results/
│       └── mock_alignments/
├── app/
│   └── streamlit_app.py
├── cli/
│   └── main.py
├── config/
│   └── default_config.yaml
├── cache/
│   └── .gitkeep
├── results/
│   └── .gitkeep
├── requirements.txt
├── pyproject.toml
├── README.md
└── DESIGN.md
```

---

## 17. Scientific Integrity Rules

### 17.1 What This System Will Never Do

| Prohibited Behavior | Reason |
|---------------------|--------|
| Return unvalidated primers | Biological safety risk |
| Use local alignment as substitute for global MSA | Scientifically incorrect |
| Silently discard BLAST hits | Loss of traceability |
| Apply heuristics without disclosure | Non-reproducible science |
| Skip secondary structure check | Off-target / failed PCR risk |
| Accept primers without specificity check | Off-target amplification risk |
| Generate results non-deterministically | Non-reproducibility |
| Truncate output for display convenience | Hidden information |

### 17.2 Disclosure Requirements

Any approximation used must appear in the output:

```
⚠ APPROXIMATION NOTICE: ΔG calculation used simplified 2-state model
  (nearest-neighbor parameters). For high-accuracy ΔG, use NUPACK or
  Vienna RNAfold. This estimate is disclosed as an approximation.
```

### 17.3 Intended Scope

This system is a **computational research aid**. It does not:
- Replace wet-lab primer validation
- Guarantee experimental PCR success
- Constitute clinical diagnostic clearance
- Account for methylation, secondary DNA structure at template level, or in-vivo conditions

---

## 18. Dependency Matrix

| Package | Version | Purpose |
|---------|---------|---------|
| Python | ≥ 3.12 | Runtime |
| Biopython | ≥ 1.83 | BLAST, FASTA parsing, SeqRecord |
| numpy | ≥ 1.26 | Numerical conservation scoring |
| pandas | ≥ 2.1 | Tabular output, CSV export |
| scipy | ≥ 1.11 | Entropy calculations |
| reportlab | ≥ 4.0 | PDF report generation |
| streamlit | ≥ 1.30 | Web interface |
| pytest | ≥ 7.4 | Test framework |
| pytest-cov | ≥ 4.1 | Coverage reporting |
| loguru | ≥ 0.7 | Structured logging |
| rich | ≥ 13.0 | Console formatting |
| pydantic | ≥ 2.5 | Data validation and type enforcement |
| click | ≥ 8.1 | CLI framework |

**External Tools (must be installed separately):**

| Tool | Version | Purpose |
|------|---------|---------|
| ClustalW2 or Clustal Omega | Latest | Global MSA (mandatory) |
| BLAST+ suite | ≥ 2.14 | Local BLASTn (optional offline mode) |

---

*End of DESIGN.md*

---

> **Scientific Disclaimer:** This system is intended for research purposes only. All primers generated must be empirically validated in a wet laboratory before use in clinical, diagnostic, or therapeutic applications. Computational predictions, however thorough, do not substitute for experimental confirmation.
