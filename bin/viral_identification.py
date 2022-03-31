import argparse
import csv

def header_file2list(header_file):
    header = []
    f = open(header_file, "r")
    content = f.read()
    header = content.splitlines()
    return header

def check_accession_species(acc, rv, hcov, hmpv, hrsv, hpiv):
    if acc in rv: return 'rv'
    elif acc in hcov: return 'hcov'
    elif acc in hmpv: return 'hmpv'
    elif acc in hrsv: return 'hrsv'
    elif acc in hpiv: return 'hpiv'
    else: return 'unidentified'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find the organism with the highest bbmap median_fold in each virus.')
    parser.add_argument('-bbmap_covstats', metavar='bbmap_covstats', required=True, type=str, help='bbmap covstats output file')
    parser.add_argument('-b', metavar='sample basename', required=True, type=str, help='sample basename')
    parser.add_argument('-rv_id', metavar='RV reference header', type=str, help='Rhinovirus reference header file.')
    parser.add_argument('-hcov_id', metavar='HCOV reference header', type=str, help='Human coronavirus reference header file.')
    parser.add_argument('-hmpv_id', metavar='HMPV reference header', type=str, help='Rhinovirus reference header file.')
    parser.add_argument('-hrsv_id', metavar='HRSV reference header', type=str, help='Rhinovirus reference header file.')
    parser.add_argument('-hpiv_id', metavar='HPIV reference header', type=str, help='Rhinovirus reference header file.')
    parser.add_argument('-m', metavar='M', type=int, default=5, help='minimum median threshold in bbmap covstats output for a reference to be considered. (Default 5)')

    args = parser.parse_args()

    rv_header = header_file2list(args.rv_id)
    hcov_header = header_file2list(args.hcov_id)
    hmpv_header = header_file2list(args.hmpv_id)
    hrsv_header = header_file2list(args.hrsv_id)
    hpiv_header = header_file2list(args.hpiv_id)

    df = []
    with open(args.bbmap_covstats,'r') as file:
        reader = csv.reader(file, delimiter="\t")
        header = next(reader)
        for i in reader:
            if int(i[9]) >= args.m: df.append(i)
    file.close()

    acc_sp_list = []
    species_list = []
    for i in range(0,len(df)):
        # add .split(" ")[0] to get just the accession if user supplies own fasta file who sequence header usually contain more than just the accessions.  
        id = df[i][0].split(" ")[0] 
        sp = check_accession_species(id, rv_header, hcov_header, hmpv_header, hrsv_header, hpiv_header)
        if not sp in species_list: 
            species_list.append(sp)
            acc_sp_list.append([id,sp])

    if len(acc_sp_list) == 0:
        output_file_name = args.b + '_failed_assembly.txt'
        output_file = open(output_file_name, 'w')
        output_file.close()
    else:
        for i,j in zip(acc_sp_list, range(1,len(acc_sp_list)+1)):
            output_file_name = args.b + '_' + i[0] + '_' + i[1]  + '_vid.txt'
            output_file = open(output_file_name, 'w+')
            output_file.write(args.b + "\t" + i[0] + "\t" + i[1] + "\n")
            output_file.close()
