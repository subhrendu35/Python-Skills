# Task 4 : Ion Statistics Calculator
# By Subhrendu Gupta (University ID : 14359440)

import argparse
import sys


def readData(file, lower, upper):
    """
    Reading the peptide mass data from task-3 output file.
    Returning a dictionary
        peps[pepid] = {
            "protein": protein,
            "pepnum": pepnum,
            "mass": mass,
            "sequence": sequence,
            
        }
        Only those peptides whose mass falls within  the m/z range given are returned
        """
    peps = {}

    # opening the file
    try:
        dataFile = open(file, "r")
    except FileNotFoundError:
        print("ERROR: Input File not found:", file)
        return peps

    # Reading the file line by line
    for line in dataFile:
        line = line.strip()

        # Skipping comments or empty lines
        if not line or line.startswith('#'):
            continue

        parts = line.split()

        # The expected components: protein, pepnum, mass, z, missedcleavage, sequence
        if len(parts) < 6:
            print("WARNING: Skipping the malformed line:", line)
            continue

        protein = parts[0]
        pepnum = parts[1]
        mass = float(parts[2])
        sequence = parts[5]

        # Using filter which is based on m/z range
        if mass < lower or mass > upper:
            continue

        # The peptide's unique identifier
        pepid = f"{protein}_{pepnum}_{sequence}"

        # Storing the peptide entry
        peps[pepid] = {
            "protein": protein,
            "pepnum": pepnum,
            "mass": mass,
            "sequence": sequence
        }

    dataFile.close()
    return peps


def main():
    # Using Argparse
    parser = argparse.ArgumentParser(
        description="Task 4: Calculating the Ion Statistics")

    # Taking  mass file as input (Task 3 output)
    parser.add_argument("fileName", help="Input peptide mass file")

    # Selecting the Range (DEFAULT = 1000-1500)
    parser.add_argument("-s",
                        "--start",
                        type=float,
                        default=1000.0,
                        help="Lower bound of m/z range")
    parser.add_argument("-e",
                        "--end",
                        type=float,
                        default=1500.0,
                        help="Upper bound of m/z range")

    # Histogram Bin Size (Mode 2)
    parser.add_argument("-b",
                        "--binsize",
                        type=float,
                        default=1.0,
                        help="Histogram Binsize (Mode 2)")

    # MODE selection
    parser.add_argument(
        "-m",
        "--mode",
        type=int,
        required=True,
        choices=[1, 2, 3, 4],
        help="1=Count, 2=Histogram, 3=Sliding Window, 4=Unique Proteins")

    # Arguments for Sliding Window (Mode 3)
    parser.add_argument("--window",
                        type=float,
                        default=200.0,
                        help="Window size (Mode 3)")
    parser.add_argument("--step",
                        type=float,
                        default=50.0,
                        help="Step size (Mode 3)")

    args = parser.parse_args()

    # Loading the Data
    file = args.fileName
    range0 = args.start
    range1 = args.end
    binsize = args.binsize

    # Loading the peptides which are filtered by m/z range
    pepmass = readData(file, range0, range1)

    # Mode 1 : Counting the Peptides in Range
    if args.mode == 1:
        print("\n== MODE 1: Peptide Count in Range ==")

        total_peptides = len(pepmass)

        # Counting the unique proteins appearing in this range
        proteins = set()
        for pid in pepmass:
            proteins.add(pepmass[pid]["protein"])

        print(f"m/z range: {range0}-{range1}")
        print(f"Total peptides: {total_peptides}")
        print(f"Total proteins represented: {len(proteins)}\n")

        sys.exit()

    # Mode 2 : Histogram Bin Size
    if args.mode == 2:
        print("\n== MODE 2: Histogram Binning ===")

        if binsize <= 0:
            print("ERROR: binsize must be > 0")
            sys.exit()

        numbins = int((range1 - range0) / binsize)

        print(f"m/z range: {range0} - {range1}")
        print(f"Binsize: {binsize}")
        print(f"Number of bins: {numbins}\n")

        # Using the Loop through bins
        for b in range(numbins):
            lower = range0 + b * binsize
            upper = lower + binsize

            count = 0
            for pid in pepmass:
                mz = pepmass[pid]["mass"]
                if lower <= mz < upper:
                    count += 1

            print(f"{lower:10.3f} - {upper:10.3f} : {count}")

        sys.exit()

    # Mode 3 : Sliding Window
    if args.mode == 3:
        print("\n== MODE 3: Analysing through Sliding Window ==")

        W = args.window
        S = args.step

        if W <= 0 or S <= 0:
            print("ERROR: window and step must be >0")
            sys.exit()

        if S >= W:
            print("ERROR: step size must be smaller than window size")
            sys.exit()

        print(f"Window size: {W}")
        print(f"Step size: {S}")
        print(f"m/z range: {range0}-{range1}\n")

        start = range0
        end = start + W

        while end <= range1:

            count = 0
            for pid in pepmass:
                mz = pepmass[pid]["mass"]
                if start <= mz < end:
                    count += 1

            print(f"{start:10.3f} - {end:10.3f} : {count}")

            start += S
            end = start + W

        sys.exit()

    # Mode 4 : Identifying the Unique Protein
    if args.mode == 4:
        print("\n== MODE 4: Unique Protein Identification ==")

        mass_to_proteins = {}

        # Creating the m/z protein-list map
        for pid in pepmass:
            mz = pepmass[pid]["mass"]
            protein = pepmass[pid]["protein"]

            if mz not in mass_to_proteins:
                mass_to_proteins[mz] = set()

            mass_to_proteins[mz].add(protein)

        # Identifying the uniquely occuring m/z values
        unique_masses = []
        for mz, prots in mass_to_proteins.items():
            if len(prots) == 1:
                unique_masses.append(mz)

        # Determining which proteins these uique masses are belonging to
        uniquely_identified = set()
        for mz in unique_masses:
            protein = list(mass_to_proteins[mz])[0]
            uniquely_identified.add(protein)

        print(f"Total peptides: {len(pepmass)}")
        print(f"Uniquely identified peptide masses: {len(unique_masses)}")
        print(
            f" The Uniquely identified proteins: {len(uniquely_identified)}\n")

        print("List of uniquely identified proteins:")
        for p in sorted(uniquely_identified):
            print(" ", p)

        sys.exit()


if __name__ == "__main__":
    main()
