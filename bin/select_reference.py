import argparse
import csv 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find the accession number associated with the highest bbmap median_fold for each virus.')
    parser.add_argument('-bbmap_covstats', metavar='bbmap_covstats', required=True, type=str, help='bbmap covstats output file')
    parser.add_argument('-b', metavar='sample basename', required=True, type=str, help='sample basename')
    parser.add_argument('-ref_header', metavar='reference header(s)', type=str, help='Reference header file.')
    parser.add_argument('-m', metavar='M', type=int, default=5, help='minimum median threshold in bbmap covstats output for a reference to be considered. (Default 5)')

    args = parser.parse_args()

    f = open(args.ref_header)
    header = []
    for i in f.read().split('\n'):
        header.append([i])

    header_dict = {}
    for i in header:
        temp = " ".join(i[0].split(' ')[1:]) 
        if not temp in header_dict:
            header_dict[temp] = [i[0].split(' ')[0]]
        else:
            header_dict[temp].append(i[0].split(' ')[0])

    map_stats = []
    with open(args.bbmap_covstats,'r') as file:
        reader = csv.reader(file, delimiter="\t")
        header = next(reader)
        for i in reader:
            if int(i[9]) >= args.m: map_stats.append(i)
        file.close()

    init_ref_header = {}
    for i in map_stats: 
    # add .split(" ")[0] to get just the accession if user supplies own fasta file who sequence header usually contain more than just the accessions.  
        sp = " ".join(i[0].split(' ')[1:])
        id = i[0].split(' ')[0]
        if not sp in init_ref_header: 
            init_ref_header[sp] = id

    if init_ref_header:
        for i in init_ref_header:
            output_file_name = args.b+'_'+str(init_ref_header[i])+'_'+'vid.txt'
            output_file = open(output_file_name, 'w+')
            output_file.write(args.b + "\t" + str(init_ref_header[i]) + "\t" + i + "\n")
            output_file.close()        
    else:
        output_file_name = args.b + '_failed_assembly.txt'
        output_file = open(output_file_name, 'w')
        map_stats_f = open(args.bbmap_covstats, 'r')
        output_file.write(map_stats_f.read())
        output_file.close()

