fasta = open("ENTER_GFACS_GENERATED_AMINO_ACID_FASTA", "r")
fasta_lines = []

for line in fasta.readlines():
    fasta_lines.append(line)

del fasta_lines[0]

number_of_sequences = 0

for element in fasta_lines:
    if element[0] == ">":
        number_of_sequences += 1

sequences = []

x = 0
counter = 0
while x < len(fasta_lines):
    if fasta_lines[x][0] != ">":
            counter = counter + 1
            sequence = fasta_lines[x].strip("\n")
            x = x + 1
            while x < len(fasta_lines):
                if fasta_lines[x][0] != ">":
                    sequence = sequence + fasta_lines[x].strip("\n")
                    x = x + 1
                else:
                    print("Just processed sequence %d out of %d (includes redundants)" % (counter,number_of_sequences))
                    sequences.append(sequence)
                    break
    else:
        x = x + 1

headers = []

sequences.append(fasta_lines[len(fasta_lines)-1])

for element in fasta_lines:
    if element[0] == ">":
        headers.append(element)

print("Removing redundants. . .\nChecking parsing")

if len(headers)==len(sequences):
    print("Numbers of headers and sequences parsed match.")
else:
    print("Press CTRL+C to cancel and check your fasta now!\Do not take next output!")

print("Please by patient while system-checking is performed. . .")

csv_fasta = []
redundant_removed_fasta_lines = []

for x in range(0, len(headers)):
    csv_fasta.append(headers[x]+","+sequences[x]+"\n")

redundant_removed_fasta_lines = []

for x in range(0, len(csv_fasta)):
    redundant_removed_fasta_lines.append(csv_fasta[x].strip("\n"))
    
redundant_removed_fasta_lines = list(set(redundant_removed_fasta_lines))

number_of_sequences = len(redundant_removed_fasta_lines)

for x in range(0,len(redundant_removed_fasta_lines)):
    redundant_removed_fasta_lines[x] = redundant_removed_fasta_lines[x].split(",")

outfile = open("C://Users/Wolf/Desktop/removed_redundants.faa", "a+")

x = 0
for x in range(0,len(redundant_removed_fasta_lines)):
    outfile.write(redundant_removed_fasta_lines[x][0]+"\n")
    outfile.write(redundant_removed_fasta_lines[x][1]+"\n")

outfile.close()

sequences2 = list(set(sequences))

redundant_seqs = []

for x in range(0, len(sequences2)):
    check = 0
    for y in range(0, len(sequences)):
        if sequences2[x] == sequences[y]:
            check = check + 1
        else:
            pass
        if check > 2:
            redundant_seqs.append(sequences2[x])
            break
        else:
            pass


redundancy = open("C://Users/Wolf/Desktop/redundancy.txt", "a+")

for element in redundant_seqs:
    for x in range(0, len(csv_fasta)):
        if csv_fasta[x].find(element+"\n") != -1:
            line = csv_fasta[x].split(",")
            redundancy.write(line[0]+"\n")
            redundancy.write(line[1]+"\n")
        else:
            pass

redundancy.close()

print("Redundancy file written")

if len(sequences) == number_of_sequences:
    print("There are %d unique protein sequences" % (number_of_sequences))
else:
    print("Multiple headers have the same sequence! Consult with Madison!")

double_partials = []
partials_5p = []
partials_3p = []
completes = []

counter_doubles = 0
counter_5p = 0
counter_3p = 0
counter_completes = 0

for sequence in sequences:
    if sequence[0] != 'M' and sequence[len(sequence)-1] != '*':
        counter_doubles = counter_doubles + 1
        print("Found double partial %d" % (counter_doubles))
        double_partials.append(sequence)
    elif sequence[0] == 'M' and sequence[len(sequence)-1] != '*':
        counter_3p = counter_3p + 1
        print("Found 3p partial %d" % (counter_3p))
        partials_3p.append(sequence)
    elif sequence[0] != 'M' and sequence[len(sequence)-1] == '*':
        counter_5p = counter_5p + 1
        print("Found 5p partial %d" % (counter_5p))
        partials_5p.append(sequence)
    elif sequence[0] == 'M' and sequence[len(sequence)-1] == '*':
        counter_completes = counter_completes + 1
        print("Found complete gene %d" % (counter_completes))
        completes.append(sequence)

print("Your breakdown is as follows. . . ")
print("Double Partials\t5p Partials\t 3p Partials\t Complete Genes\t Total")
print("%d\t%d\t%d\t%d\t%d" % (counter_doubles, counter_5p, counter_3p, counter_completes, counter_doubles+counter_5p+counter_3p+counter_completes))
print("Writing breakdown to files. . .")

three_primes = open("C://Users/Wolf/Desktop/three_primes.txt", "a+")
five_primes = open("C://Users/Wolf/Desktop/five_primes.txt", "a+")
doubles = open("C://Users/Wolf/Desktop/doubleprimes.txt", "a+")
completes_dest = open("C://Users/Wolf/Desktop/completes.txt", "a+")

for element in double_partials:
     for x in range(0, len(csv_fasta)):
        if csv_fasta[x].find(element+"\n") != -1:
            line = csv_fasta[x].split(",")
            doubles.write(line[0]+"\n")
            doubles.write(line[1]+"\n")
        else:
            pass    
doubles.close()
print("Doubles written")

for element in partials_5p:
    for x in range(0, len(csv_fasta)):
        if csv_fasta[x].find(element+"\n") != -1:
            line = csv_fasta[x].split(",")
            five_primes.write(line[0]+"\n")
            five_primes.write(line[1]+"\n")
        else:
            pass  

five_primes.close()
print("5p written")

for element in partials_3p:
    for x in range(0, len(csv_fasta)):
        if csv_fasta[x].find(element+"\n") != -1:
            line = csv_fasta[x].split(",")
            three_primes.write(line[0]+"\n")
            three_primes.write(line[1]+"\n")
        else:
            pass

three_primes.close()
print("3P written")

for element in list(set(completes)):
    for x in range(0, len(csv_fasta)):
        if csv_fasta[x].find(element+"\n") != -1:
            line = csv_fasta[x].split(",")
            completes_dest.write(line[0]+"\n")
            completes_dest.write(line[1]+"\n")
        else:
            pass

completes_dest.close()
print("Completes written")

print("Goodbye")
