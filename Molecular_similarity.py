import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from pathlib import Path

def calculate_molecular_similarity(reference_smiles, db_input, output_results):

    # 1. Initialize reference molecule
    ref_mol = Chem.MolFromSmiles(reference_smiles)

    if ref_mol is None:
        raise ValueError(
            f"Failed to decode reference SMILES: {reference_smiles}"
        )

    # Generate the reference Morgan fingerprint
    # Radius = 2, 2048 bits, pharmacophoric features enabled
    ref_fp = AllChem.GetMorganFingerprintAsBitVect(
        ref_mol,
        2,
        nBits=2048,
        useFeatures=True
    )

    # 2. Read input database
    input_file = Path(db_input)
    output_file = Path(output_results)

    try:
        # Expected input format:
        # Column 1: CAS Registry Number
        # Column 2: SMILES
        # Tab-delimited, without a header
        df_library = pd.read_csv(
            input_file,
            sep="\t",
            names=["CAS", "SMILES"],
            header=None,
            dtype=str
        )

    except Exception as e:
        print(f"Error accessing input file: {e}")
        return

    processed_data = []

    print(
        f"Executing similarity screening on "
        f"{len(df_library)} compounds..."
    )

    # 3. Calculate molecular similarity
    for _, row in df_library.iterrows():

        cas_id = row["CAS"]

        target_smi = (
            str(row["SMILES"]).strip()
            if pd.notnull(row["SMILES"])
            else ""
        )

        # Skip empty SMILES entries
        if not target_smi:
            continue

        # Convert SMILES to an RDKit molecule
        current_mol = Chem.MolFromSmiles(target_smi)

        # Skip invalid SMILES
        if current_mol is None:
            continue

        # Generate Morgan fingerprint
        current_fp = AllChem.GetMorganFingerprintAsBitVect(
            current_mol,
            2,
            nBits=2048,
            useFeatures=True
        )

        # Calculate Tanimoto similarity
        tanimoto_score = DataStructs.TanimotoSimilarity(
            ref_fp,
            current_fp
        )

        # Store results for all valid compounds
        processed_data.append(
            {
                "CAS_Registry": cas_id,
                "SMILES_String": target_smi,
                "Tanimoto_Similarity": round(tanimoto_score, 4)
            }
        )

    # 4. Export results
    if processed_data:

        results_df = pd.DataFrame(processed_data)

        results_df.to_csv(
            output_file,
            sep="\t",
            index=False
        )

        print(
            f"Analysis complete. Results stored in: "
            f"{output_file}"
        )

    else:
        print("No valid structures found to analyze.")


if __name__ == "__main__":

    # Reference parent neonicotinoid
    NEO_SMILES = (
        "C1CN(C(=N[N+](=O)[O-])N1)"
        "CC2=CN=C(C=C2)Cl"
    )

    # File configuration
    INPUT_DB = "Existing_Chemical_Database.txt"
    OUTPUT_NAME = "Similarity_Result.txt"

    # Run structural similarity analysis
    calculate_molecular_similarity(
        NEO_SMILES,
        INPUT_DB,
        OUTPUT_NAME
    )