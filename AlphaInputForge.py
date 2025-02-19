"""
AlphaFold 3 JSON Generator with MMseqs2 MSA (Single Database Support)

This script:
1. Searches a single MMseqs2 database to generate an unpaired MSA.
2. Uses 'align -a' to ensure backtrace CIGAR is included for AlphaFold 3.
3. Uses '--format-mode 3' to produce a correct .a3m file.
4. Loads ligands from corresponding .tsv files (or Uniform.tsv as fallback) and JSON-escapes the SMILES using jq.
5. Generates JSON input files for AlphaFold 3 with relative paths for the .a3m file.
6. Orders the JSON so that all protein entries appear first and then any ligand entries.
7. If no ligand data exists for a protein, no ligand entry is added.

USAGE:
python generate_alphafold3_json.py --input_path path_to_input_folder \
--output_path path_to_output_folder --cpu 16 --mmseqs_bin path_to/mmseqs \
--mmseqs_DB path_to/uniref90 --verbose
"""

import os
import json
import argparse
import subprocess
import datetime
from tqdm import tqdm
from Bio import SeqIO
import shutil

class Logger:
    """Handles logging of messages, warnings, and errors with verbosity control."""
    def __init__(self, log_file, verbose=False):
        self.log_file = log_file
        self.verbose = verbose

    def log(self, message, is_error=False):
        """Logs a message with a timestamp."""
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_entry = f"[{timestamp}] {message}"
        if self.verbose or is_error:
            print(log_entry)
        with open(self.log_file, "a") as log:
            log.write(log_entry + "\n")
        if is_error:
            print("⚠️ ERROR:", message)

class MSAGenerator:
    """Handles MSA generation using a single MMseqs2 database."""
    def __init__(self, num_threads, output_dir, mmseqs_path, db_path, logger):
        self.num_threads = num_threads
        self.output_dir = output_dir
        self.mmseqs_path = mmseqs_path
        self.db_path = db_path
        self.logger = logger
        self.tmp_path = os.path.join(output_dir, "tmp")
        os.makedirs(self.tmp_path, exist_ok=True)

    def generate_msa(self, protein_id, fasta_file):
        """Runs MMseqs2 to generate an unpaired MSA."""
        msa_file = os.path.join(self.output_dir, f"{protein_id}_unpaired.a3m")
        db_local = os.path.join(self.output_dir, f"{protein_id}_db")
        result_path = os.path.join(self.output_dir, f"{protein_id}_result")
        align_path = os.path.join(self.output_dir, f"{protein_id}_aligned")

        try:
            # Create MMseqs2 database from the FASTA file
            subprocess.run([self.mmseqs_path, "createdb", fasta_file, db_local], check=True)

            # Run MMseqs2 search
            subprocess.run([
                self.mmseqs_path, "search", db_local, self.db_path, result_path,
                self.tmp_path, "--threads", str(self.num_threads), "--max-seqs", "1000"
            ], check=True)

            # Run MMseqs2 align with the -a flag to include backtrace CIGAR
            subprocess.run([
                self.mmseqs_path, "align", db_local, self.db_path, result_path,
                align_path, "--threads", str(self.num_threads), "-a"
            ], check=True)

            # Convert alignment to A3M format
            subprocess.run([
                self.mmseqs_path, "convertalis", db_local, self.db_path, align_path,
                msa_file, "--format-mode", "3"
            ], check=True)

            return msa_file if os.path.exists(msa_file) else ""

        except subprocess.CalledProcessError:
            self.logger.log(f"❌ MMseqs2 failed for {protein_id}. Skipping MSA.", is_error=True)
            return ""

class JSONGenerator:
    """Generates JSON input files for AlphaFold 3 from FASTA and TSV files."""
    def __init__(self, input_folder, output_folder, num_threads, mmseqs_path, db_path, verbose):
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.num_threads = num_threads
        self.logger = Logger(os.path.join(output_folder, "process.log"), verbose)
        self.msa_generator = MSAGenerator(num_threads, output_folder, mmseqs_path, db_path, self.logger)
        os.makedirs(self.output_folder, exist_ok=True)

    def load_ligand_data(self, tsv_file):
        """Loads ligand data from a TSV file into a dictionary, omitting empty lines."""
        ligand_dict = {}
        if os.path.exists(tsv_file):
            with open(tsv_file, "r") as file:
                for line in file:
                    if not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 3:
                        self.logger.log(f"⚠️ Skipping invalid line in {tsv_file}: {line.strip()}", is_error=True)
                        continue
                    protein_id, ligand_id, smiles = parts[:3]
                    try:
                        escaped_smiles = subprocess.check_output(
                            ["jq", "-R", "."],
                            input=smiles,
                            text=True
                        ).strip()
                        if escaped_smiles.startswith('"') and escaped_smiles.endswith('"'):
                            escaped_smiles = escaped_smiles[1:-1]
                    except subprocess.CalledProcessError as e:
                        self.logger.log(f"⚠️ Error JSON-escaping SMILES: {e}", is_error=True)
                        escaped_smiles = smiles
                    if protein_id not in ligand_dict:
                        ligand_dict[protein_id] = []
                    ligand_dict[protein_id].append({"ligand": {"id": [ligand_id], "smiles": escaped_smiles}})
        return ligand_dict

    def process_fasta_file(self, fasta_file):
        """Processes a FASTA file, generates MSA, and creates a JSON file with ligands.
        
        All protein entries are added first, followed by all ligand entries (if any).
        """
        fasta_path = os.path.join(self.input_folder, fasta_file)
        json_file = os.path.join(self.output_folder, f"{fasta_file}.json")
        proteins = []
        ligands = []

        # Determine the corresponding ligand TSV file; fallback to Uniform.tsv if not found
        tsv_file = os.path.join(self.input_folder, fasta_file.replace(".fasta", ".tsv"))
        if not os.path.exists(tsv_file):
            tsv_file = os.path.join(self.input_folder, "Uniform.tsv")
        ligand_dict = self.load_ligand_data(tsv_file) if os.path.exists(tsv_file) else {}

        self.logger.log(f"🚀 Processing {fasta_file}")
        fasta_sequences = list(SeqIO.parse(fasta_path, "fasta"))

        for record in tqdm(fasta_sequences, desc=f"Processing {fasta_file}", unit="seq"):
            protein_id = record.id
            sequence = str(record.seq)
            unpaired_msa_abs = self.msa_generator.generate_msa(protein_id, fasta_path)
            unpaired_msa_rel = os.path.relpath(unpaired_msa_abs, start=self.output_folder) if unpaired_msa_abs else ""

            protein_entry = {"id": [protein_id], "sequence": sequence}
            if unpaired_msa_rel:
                protein_entry["unpairedMsaPath"] = unpaired_msa_rel
            proteins.append({"protein": protein_entry})

            # Only add ligand entries if available; omit otherwise.
            if protein_id in ligand_dict:
                ligands.extend(ligand_dict[protein_id])

        # Final JSON: All protein entries first, then all ligand entries.
        sequences = proteins + ligands

        with open(json_file, "w") as jf:
            json.dump({
                "name": fasta_file,
                "sequences": sequences,
                "dialect": "alphafold3",
                "version": 2
            }, jf, indent=4)

        self.logger.log(f"✅ Completed {fasta_file} -> {json_file}")

    def run(self):
        """Processes all FASTA files in the input directory."""
        fasta_files = [f for f in os.listdir(self.input_folder) if f.endswith(".fasta")]
        if not fasta_files:
            self.logger.log("⚠️ No FASTA files found!", is_error=True)
            return
        for fasta_file in tqdm(fasta_files, desc="Processing FASTA files", unit="file"):
            self.process_fasta_file(fasta_file)
        self.logger.log("✅ All FASTA files processed successfully!")

def prepare_output_folder(output_path):
    """Checks if the output folder exists and asks whether to delete it or create a new one."""
    if os.path.exists(output_path):
        response = input(f"Output folder '{output_path}' already exists. Delete it and create a new one? (y/n): ")
        if response.lower() == 'y':
            try:
                shutil.rmtree(output_path)
                os.makedirs(output_path)
                print(f"Deleted and recreated output folder: {output_path}")
            except Exception as e:
                print(f"❌ Failed to delete output folder: {e}")
                exit(1)
        else:
            # Create a new output folder by appending a timestamp
            new_output = output_path.rstrip("/\\") + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            os.makedirs(new_output)
            print(f"Using new output folder: {new_output}")
            return new_output
    return output_path

if __name__ == "__main__":
    import shutil
    parser = argparse.ArgumentParser(description="Generate JSON inputs for AlphaFold 3 using MMseqs2.")
    parser.add_argument("--input_path", type=str, required=True, help="Path to input folder containing FASTA & TSV files")
    parser.add_argument("--output_path", type=str, required=True, help="Path to output folder")
    parser.add_argument("--cpu", type=int, required=True, help="Number of CPU threads")
    parser.add_argument("--mmseqs_bin", type=str, required=True, help="Path to MMseqs2 binary")
    parser.add_argument("--mmseqs_DB", type=str, required=True, help="Path to the single MMseqs2 indexed database")
    parser.add_argument("--verbose", action="store_true", help="Enable detailed logging")
    args = parser.parse_args()

    output_folder = prepare_output_folder(args.output_path)
    JSONGenerator(args.input_path, output_folder, args.cpu, args.mmseqs_bin, args.mmseqs_DB, args.verbose).run()
