"""
Test script for CPTAC Browser data loading functions
Run this before deploying to verify all data loaders work correctly
"""

import cptac
import pandas as pd
import numpy as np
from typing import Optional, Tuple, List

def process_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Handle multi-index and transpose dataframe"""
    print(f"  Input type: {type(df)}, shape: {df.shape if hasattr(df, 'shape') else 'N/A'}")

    if isinstance(df.columns, pd.MultiIndex):
        print(f"  Multi-index detected with {df.columns.nlevels} levels")
        # Deduplicate multi-index columns by taking mean
        df = df.groupby(level=list(range(df.columns.nlevels)), axis=1).mean()
        print(f"  After groupby type: {type(df)}")

    # Ensure we have a proper DataFrame
    if not isinstance(df, pd.DataFrame):
        print(f"  Converting {type(df)} to DataFrame")
        df = pd.DataFrame(df)

    result = df.T
    print(f"  Output type: {type(result)}, shape: {result.shape}")
    return result

def load_phosphoproteomics(data) -> Optional[pd.DataFrame]:
    """Load phosphoproteomics data with fallback to multiple sources"""
    sources = ['harmonized', 'broad', 'mssm', 'washu', None]

    for source in sources:
        try:
            print(f"    Trying source: {source}")
            if source:
                phospho_raw = data.get_phosphoproteomics(source=source)
            else:
                phospho_raw = data.get_phosphoproteomics()

            if phospho_raw is not None and not phospho_raw.empty:
                print(f"    ✓ Success with source: {source}")
                return phospho_raw
        except Exception as e:
            print(f"    ✗ Failed with source {source}: {e}")
            continue

    return None

def test_cancer_type(cancer_code: str, cancer_name: str):
    """Test data loading for a specific cancer type"""
    print(f"\n{'='*60}")
    print(f"Testing: {cancer_name} ({cancer_code})")
    print(f"{'='*60}")

    try:
        # Load cancer data
        print(f"\n1. Loading cancer data...")
        cancer_map = {
            "brca": cptac.Brca,
            "coad": cptac.Coad,
            "hnscc": cptac.Hnscc,
            "luad": cptac.Luad,
            "ov": cptac.Ov,
            "ccrcc": cptac.Ccrcc,
            "gbm": cptac.Gbm,
            "lscc": cptac.Lscc,
            "pdac": cptac.Pdac
        }
        data = cancer_map[cancer_code]()
        print(f"  ✓ Cancer data loaded")

        # Test proteomics
        print(f"\n2. Testing proteomics data...")
        try:
            proteomics_raw = data.get_proteomics()
            print(f"  Raw proteomics shape: {proteomics_raw.shape}")
            proteomics = process_dataframe(proteomics_raw)
            print(f"  ✓ Proteomics processed successfully")
            print(f"    Genes: {len(proteomics.index)}")
            print(f"    Samples: {len(proteomics.columns)}")
            print(f"    First 3 genes: {list(proteomics.index[:3])}")
        except Exception as e:
            print(f"  ✗ Proteomics failed: {e}")
            import traceback
            traceback.print_exc()

        # Test phosphoproteomics
        print(f"\n3. Testing phosphoproteomics data...")
        try:
            phospho_raw = load_phosphoproteomics(data)
            if phospho_raw is not None:
                print(f"  Raw phospho shape: {phospho_raw.shape}")
                phospho = process_dataframe(phospho_raw)
                print(f"  ✓ Phospho processed successfully")
                print(f"    Phosphosites: {len(phospho.index)}")
                print(f"    Samples: {len(phospho.columns)}")
                print(f"    First 3 sites: {list(phospho.index[:3])}")
            else:
                print(f"  ⚠ No phosphoproteomics data available")
        except Exception as e:
            print(f"  ✗ Phospho failed: {e}")
            import traceback
            traceback.print_exc()

        print(f"\n✓ {cancer_name} tests completed")
        return True

    except Exception as e:
        print(f"\n✗ {cancer_name} failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run tests on all cancer types"""
    print("CPTAC Browser Data Loading Test Suite")
    print("=" * 60)

    cancer_types = {
        "brca": "Breast Cancer",
        "coad": "Colon Adenocarcinoma",
        "ov": "Ovarian Cancer",
        # Add more as needed
    }

    results = {}
    for code, name in cancer_types.items():
        results[name] = test_cancer_type(code, name)

    # Summary
    print(f"\n{'='*60}")
    print("TEST SUMMARY")
    print(f"{'='*60}")
    for name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {name}")

    total = len(results)
    passed = sum(results.values())
    print(f"\nTotal: {passed}/{total} passed")

    return all(results.values())

if __name__ == "__main__":
    import sys
    success = main()
    sys.exit(0 if success else 1)
