
---

# AlphaInputForge

**AlphaInputForge** is a command‐line tool that streamlines the preparation of input files for AlphaFold 3. It automatically generates multiple sequence alignments (MSAs) using MMseqs2 (with alignment backtraces enabled), processes ligand data from TSV files, and produces well‐formatted JSON files that conform to AlphaFold 3’s input specifications.

The tool supports:
- **Automated MSA Generation:**  
  Uses MMseqs2 to generate unpaired MSAs in A3M format, ensuring that alignment backtrace information is included (using the `-a` flag with MMseqs2 `align`).

- **Ligand Data Integration:**  
  Loads ligand data from TSV files (named to match the FASTA file) or from a fallback `Uniform.tsv` file.  
  Each TSV file should have three tab‐separated columns:  
  `Protein_ID`, `Ligand_ID`, and `SMILES` (the SMILES string is JSON-escaped using the jq command-line tool).  
  The JSON output includes all protein entries first, followed by ligand entries (if any).  
  If no ligand data exists for a protein, no ligand entry is added.

- **Output Folder Management:**  
  If the specified output folder already exists, the script will ask whether to delete it and create a new one or to create a new folder (with a timestamp appended).

- **Detailed Logging and Verbose Progress:**  
  Uses a Logger class to write timestamped messages and errors to a log file and (optionally) the console, with progress bars provided by `tqdm`.

## Features

- **Automated MSA Generation**  
  Generates MSAs for each protein in the FASTA files by converting them to an MMseqs2 database and then running MMseqs2’s search, align (with `-a`), and convertalis steps.  
  The output MSA (in A3M format) is stored with a relative path in the JSON file.

- **Ligand Data Processing**  
  Loads ligand data from TSV files where each row contains the `Protein_ID` (which must exactly match the FASTA header), `Ligand_ID`, and the SMILES string.  
  The SMILES strings are JSON-escaped using the jq command-line tool so that they can be safely embedded in JSON.

- **Flexible Input/Output**  
  Processes multiple FASTA files from a given input folder and creates a JSON file for each FASTA file.  
  Each JSON file contains all protein entries (with MSA paths) first, then any ligand entries.

- **Robust Error Handling & Logging**  
  Verbose logging with timestamps is provided, and errors are reported if, for example, MMseqs2 fails or the input files are invalid.

- **Output Folder Management**  
  If the output folder already exists, the script prompts you to either delete it or create a new folder with a timestamp.

## Installation

### System Requirements
- **Python 3.6+**
- **jq:** A lightweight and flexible command-line JSON processor.

### Installing Python Dependencies
Install the required Python packages using pip:
```bash
pip install biopython tqdm
```

### Installing jq
- **Ubuntu/Debian:**
  ```bash
  sudo apt-get update
  sudo apt-get install jq
  ```
- **Fedora/CentOS/RHEL:**
  ```bash
  sudo dnf install jq
  ```
- **macOS (using Homebrew):**
  ```bash
  brew install jq
  ```

### Installing MMseqs2
Download and install MMseqs2 from the [official website](https://mmseqs.com/). For example:
Install the right package for your system otherwise, mmseqs might not work
```bash
wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
tar xvf mmseqs-linux-avx2.tar.gz
sudo cp mmseqs/bin/mmseqs /usr/local/bin/
```
Test the installation with:
```bash
mmseqs --version
```

### Pre-indexing the Database
AlphaInputForge requires a pre-indexed MMseqs2 database. To create one:
```bash
mmseqs createdb uniref90.fasta uniref90_db
mmseqs createindex uniref90_db tmp
```
Replace `uniref90.fasta` with your database FASTA file. The resulting database (`uniref90_db`) is what you pass to the script.

## Input File Format

### FASTA Files
- **Filename Example:** `cyp1.fasta`
- **Header Example:**
  ```
  >sp|P05108|CP11A_HUMAN Cholesterol side-chain cleavage enzyme, mitochondrial OS=Homo sapiens OX=9606 GN=CYP11A1 PE=1 SV=2
  MLAKGLPPRSVLVKGCQTFLSAPREGLGRLRVPTGEGAGISTRSPRPFNEIPSPGDNGWL...
  ```
- **Notes:**  
  The FASTA header should include the protein ID (e.g., `sp|P05108|CP11A_HUMAN`).  
  The sequences should contain only valid amino acid characters.

### TSV Files
- **Filename Convention:**  
  The TSV file should have the same name as the FASTA file, but with the `.tsv` extension (e.g., `cyp1.tsv`).  
  Alternatively, if no such file exists, a file named `Uniform.tsv` will be used as a fallback.
- **File Format:**  
  Tab-separated with at least three columns:  
  - **Protein_ID:** Must match the protein ID from the FASTA header (e.g., `sp|P05108|CP11A_HUMAN`)
  - **Ligand_ID:** e.g., `Heme_B`
  - **SMILES:** e.g., `CC1=C(C2=CC3=C(C(=C([N-]3)C=C4C(=C(C(=N4)C=C5C(=C(C(=N5)C=C1[N-]2)C)C=C)C)C=C)C)CCC(=O)[O-])CCC(=O)[O-].[Fe]`
- **Example:**
  ```
  sp|P05108|CP11A_HUMAN	Heme_B	CC1=C(C2=CC3=C(C(=C([N-]3)C=C4C(=C(C(=N4)C=C5C(=C(C(=N5)C=C1[N-]2)C)C=C)C)C=C)C)CCC(=O)[O-])CCC(=O)[O-].[Fe]
  sp|P05108|CP11A_HUMAN	22-Hydroxy-Cholesterol	C[C@@H]([C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(CCC(C)C)O
  ```

## Usage

Run the script using the following command:
```bash
python generate_alphafold3_json.py --input_path ./input --output_path ./output \
  --cpu 16 --mmseqs_bin /usr/local/bin/mmseqs --mmseqs_DB ./db/uniref90 --verbose
```

### Command-Line Arguments:
- **`--input_path`**: Path to the folder containing FASTA (and optionally TSV) files.
- **`--output_path`**: Path to the folder where the JSON and MSA output will be saved.
- **`--cpu`**: Number of CPU threads to use for MMseqs2 processing.
- **`--mmseqs_bin`**: Path to the MMseqs2 binary.
- **`--mmseqs_DB`**: Path to the pre-indexed MMseqs2 database.
- **`--verbose`**: Enable detailed logging and progress messages.

## Output
- **MSA Files:**  
  Unpaired MSA files in A3M format (e.g., `cyp1_unpaired.a3m`) are generated in the output folder.  
  The JSON file stores the relative path to the MSA file.

- **JSON Files:**  
  A JSON file is created for each FASTA file. It includes:
  - **Protein Entries:**  
    Contain the protein ID, sequence, and relative path to the unpaired MSA.
  - **Ligand Entries:**  
    Each ligand entry (if present) is appended after the protein entries.
  - **Example JSON Structure:**
  ```json
  {
    "name": "cyp1.fasta",
    "sequences": [
      {
        "protein": {
          "id": ["sp|P05108|CP11A_HUMAN"],
          "sequence": "MLAKGLPPRSVLVKGCQTFLSAPREGLGRLRVPTGEGAGISTRSPRPFNEIPSPGDNGWL...",
          "unpairedMsaPath": "cyp1_unpaired.a3m"
        }
      },
      {
        "ligand": {
          "id": ["Heme_B"],
          "smiles": "CC1=C(C2=CC3=C(C(=C([N-]3)C=C4C(=C(C(=N4)C=C5C(=C(C(=N5)C=C1[N-]2)C)C=C)C)C=C)C)CCC(=O)[O-])CCC(=O)[O-].[Fe]"
        }
      },
      {
        "ligand": {
          "id": ["22-Hydroxy-Cholesterol"],
          "smiles": "C[C@@H]([C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(CCC(C)C)O"
        }
      }
    ],
    "dialect": "alphafold3",
    "version": 2
  }
  ```

## Output Folder Management
If the output folder already exists, the script will prompt:
```
Output folder 'OUTPUT_FOLDER' already exists. Delete it and create a new one? (y/n):
```
- If **`y`** is entered, the folder is deleted and recreated.
- If **`n`** is entered, a new folder with a timestamp appended to its name is created, and that folder is used as the output folder.

## Contributing
Contributions and improvements are welcome! Please open issues or submit pull requests for enhancements or bug fixes.

## License
This project is licensed under the MIT License.

---

### **Summary**
AlphaInputForge automates the generation of input files for AlphaFold 3 by:
- Generating MSAs using MMseqs2 from protein FASTA files.
- Integrating ligand data from TSV files (with JSON-escaped SMILES strings via jq).
- Creating well-structured JSON input files that list protein entries first, followed by ligand entries.
- Providing detailed logging and robust error handling.
- Offering flexible output folder management.

Feel free to modify this README to fit your needs and add any additional instructions or details!

---

This README can now be added to your GitHub repository for **AlphaInputForge**. Let me know if you need further modifications or additional sections!
