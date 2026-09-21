import pandas as pd
from rdkit import Chem
from pathlib import Path

# Predefined NEO-related structural patterns
NEO_STRUCTURAL_PATTERNS = [
    "ClC1=NC=C(C)C=C1",
    "ClC1=NC=C(C)S1",
    "CN/C(NC)=N/[N+]([O-])=O",
    "CN(CCC1)C1C2=CC=CN=C2",
    "FC(F)(F)C1=CC=NC=C1",
    "C1=CC=NC=C1"
]

def process_chemical_library(input_path: str, output_path: str):

    # 1. Set input and output paths
    source_file = Path(input_path)
    result_file = Path(output_path)

    # 2. Initialize RDKit query molecules
    patterns = []

    for smi in NEO_STRUCTURAL_PATTERNS:

        pattern = Chem.MolFromSmiles(smi)

        if pattern is None:
            raise ValueError(
                f"Invalid predefined structural pattern: {smi}"
            )

        patterns.append(pattern)

    # 3. Load chemical database
    try:
        raw_df = pd.read_csv(
            source_file,
            usecols=["CAS", "smiles"],
            dtype=str
        )

    except Exception as e:
        print(f"Error accessing input file: {e}")
        return

    matched_records = []

    print(f"Analyzing {len(raw_df)} entries...")

    # 4. Perform substructure screening
    for _, row in raw_df.iterrows():

        cas_id = row["CAS"]
        smi_content = row["smiles"]

        # Skip empty or invalid SMILES entries
        if not isinstance(smi_content, str) or not smi_content.strip():
            continue

        smi_content = smi_content.strip()

        # Convert SMILES to RDKit molecule
        target_mol = Chem.MolFromSmiles(smi_content)

        if target_mol is None:
            continue

        # Determine whether the molecule contains at least one
        # predefined NEO-related structural pattern
        match_found = any(
            target_mol.HasSubstructMatch(pattern)
            for pattern in patterns
        )

        # Store matched compounds
        if match_found:

            matched_records.append(
                {
                    "CAS": cas_id,
                    "smiles": smi_content
                }
            )

    # 5. Export matched compounds
    if matched_records:

        final_output = pd.DataFrame(matched_records)

        final_output.to_csv(
            result_file,
            index=False
        )

        print(
            f"Exported {len(final_output)} matches to "
            f"{result_file}"
        )

    else:
        print("No valid matches identified.")


if __name__ == "__main__":

    # File configuration
    INPUT_DB = "Existing_Chemical_Database.csv"
    OUTPUT_NAME = "Extraction_Result.csv"

    # Run substructure-based screening

    process_chemical_library(
        INPUT_DB,
        OUTPUT_NAME
    )