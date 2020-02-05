# imported packages/modules
from Bio import SeqIO
import collections
import re
import pprint
import time
import numpy as np
import math
from collections import defaultdict
import os

# A class to help put genome and gene annotations into data structures
# Contains methods that help find various data in the genome such as codon frequency, nucleotide frequency etc etc

class Genome_Data(object):

    # Class Constructor with given fasta file(geonome) and gff3(gene annotation),
    # saves into data structures(dictionary for genome and list for annotation)

    def __init__(self, fasta, gff3):
        self.seq_dict = {
            rec.id: rec.seq for rec in SeqIO.parse(fasta, "fasta")}
        self.seq_dict = collections.OrderedDict(sorted(self.seq_dict.items()))

        self.gene_annotation = []
        with open(gff3) as fp:
            line = fp.readline()
            cnt = 1
            while line:
                if len(line.split()) > 7:
                    if line.split()[2] == "CDS":
                        if(line.split()[6] == "+"):
                            self.gene_annotation.append(line.split())
                line = fp.readline()
                cnt += 1

    # divides and saves intergenetic regions and genic regions in the genome from given gene annotation

    def Integernetic_and_Genetic_Region(self):

        self.intergenic_region = []
        self.genic_region = []

        index = 0
        start = 0

        DNA = self.gene_annotation[index][0]

        while index < len(self.gene_annotation):

            if DNA != self.gene_annotation[index][0]:
                genome = genome = (self.seq_dict.get(DNA, ""))
                intergenic = genome[start:len(genome)]
                self.intergenic_region.append(intergenic)
                start = 0
                DNA = self.gene_annotation[index][0]

            gene_start = int(self.gene_annotation[index][3])
            gene_stop = int(self.gene_annotation[index][4])

            genome = (self.seq_dict.get(DNA, ""))

            intergenic = genome[start:gene_start-1]
            genic = genome[gene_start-1:gene_stop]

            self.intergenic_region.append(intergenic)
            self.genic_region.append(genic)

            start = gene_stop
            index = index+1

            if(index == len(self.gene_annotation)):
                intergenic = genome[start:len(genome)]
                self.intergenic_region.append(intergenic)

    # returns list of intergenic regions

    def getIntegernetic(self):
        return self.intergenic_region

    # returns list of genic regions

    def getGenetic(self):
        return self.genic_region

    # method that prints average length of intergenic region and genic region

    def printAverageLengths(self):
        self.intergenic_length = (
            sum(map(len, self.intergenic_region))/int(len(self.intergenic_region)))
        self.genic_length = (
            sum(map(len, self.genic_region))/int(len(self.genic_region)))

        print("Average Intergenic Length: {} ".format(
            sum(map(len, self.intergenic_region))/int(len(self.intergenic_region))))
        print("Average Genic Length: {} ".format(
            sum(map(len, self.genic_region))/int(len(self.genic_region))))

    # method that saves frequency tables( in dictionary form) of intergenic nucleotide frequencies,
    # genic start/middle/stop codon frequency, and genic codon frequency

    def updateFrequencyTables(self, intergenic_region, genic_region):

        A = 0
        C = 0
        G = 0
        T = 0
        total = 0

        self.start_codons = []
        self.middle_codons = []
        self.stop_codons = []

        self.nucleotide_frequency = dict()
        self.frequency_table = dict()
        self.startCodon_frequency = dict()
        self.middleCodon_frequency = dict()
        self.stopCodon_frequency = dict()

        for x in intergenic_region:
            A = A+x.count('A')
            C = C+x.count('C')
            G = G+x.count('G')
            T = T+x.count('T')
            total = total+len(x)

        self.nucleotide_frequency.update({'A': float(A)/total})
        self.nucleotide_frequency.update({'C': float(C)/total})
        self.nucleotide_frequency.update({'G': float(G)/total})
        self.nucleotide_frequency.update({'T': float(T)/total})

        total_codon = 0
        total_start_codon = 0
        total_middle_codon = 0
        total_stop_codon = 0

        for gene in genic_region:
            start = 0
            stop = len(re.findall('...', str(gene)))

            for codon in (re.findall('...', str(gene))):

                if start == 0:
                    if codon in self.startCodon_frequency:

                        self.startCodon_frequency[codon] = self.startCodon_frequency[codon]+1
                        total_start_codon = total_start_codon+1

                    else:
                        self.startCodon_frequency.update({codon: 1})
                        total_start_codon = total_start_codon+1

                elif start == stop-1:

                    if codon in self.stopCodon_frequency:

                        self.stopCodon_frequency[codon] = self.stopCodon_frequency[codon]+1
                        total_stop_codon = total_stop_codon+1

                    else:
                        self.stopCodon_frequency.update({codon: 1})
                        total_stop_codon = total_stop_codon+1
                else:
                    if codon in self.middleCodon_frequency:

                        self.middleCodon_frequency[codon] = self.middleCodon_frequency[codon]+1
                        total_middle_codon = total_middle_codon+1

                    else:
                        self.middleCodon_frequency.update({codon: 1})
                        total_middle_codon = total_middle_codon+1

                if codon in self.frequency_table:

                    self.frequency_table[codon] = self.frequency_table[codon]+1
                    total_codon = total_codon+1
                else:
                    self.frequency_table.update({codon: 1})
                    total_codon+total_codon+1

                start = start+1

        for keys in self.frequency_table.keys():
            self.frequency_table[keys] = float(
                self.frequency_table[keys])/total_codon

        for keys in self.startCodon_frequency.keys():
            self.startCodon_frequency[keys] = float(
                self.startCodon_frequency[keys])/total_start_codon

        for keys in self.stopCodon_frequency.keys():
            self.stopCodon_frequency[keys] = float(
                self.stopCodon_frequency[keys])/total_stop_codon

        for keys in self.middleCodon_frequency.keys():
            self.middleCodon_frequency[keys] = float(
                self.middleCodon_frequency[keys])/total_middle_codon

    # method that prints start codon frequency table (Genic)

    def printStartCodonFrequency(self):

        pprint.pprint(self.startCodon_frequency)

     # method that prints stop codon frequency table (Genic)

    def printStopCodonFrequency(self):

        pprint.pprint(self.stopCodon_frequency)

     # method that prints middle codon frequency table (Genic)

    def printMiddleCodonFrequency(self):

        pprint.pprint(self.middleCodon_frequency)

    # method that prints nucleotide frequency table (Genic)

    def printNucleotide_FrequencyTable(self):

        A = 0
        C = 0
        G = 0
        T = 0
        total = 0

        for x in self.intergenic_region:
            A = A+x.count('A')
            C = C+x.count('C')
            G = G+x.count('G')
            T = T+x.count('T')
            total = total+len(x)

        print("A :%f , C: %f, G: %f, T: %f " %
              (float(A)/total, float(C)/total, float(G)/total, float(T)/total))

    # method that prints codon frequency table for all codon (Genic)

    def printCodon_FrequencyTable(self):
        pprint.pprint(self.frequency_table)

    # method that compares a set of predicted genes with a set of gene annotations,(predicted genes are found with Hidden Markov Model)
    # it prints out the number of genes that have same start index, stop index, and both same and start index

    def compareWithPrediction(self, gff3):

        self.predictions = []

        with open(gff3) as fp:
            line = fp.readline()
            while line:
                self.predictions.append(line.split())
                line = fp.readline()

            end = 0
            start = 0
            both = 0

            length = len(self.gene_annotation)

            for keys in self.seq_dict:
                key = keys
                #print("for key {}".format(key))
                for x in range(length):
                    if self.gene_annotation[x][0] == key:
                        # print("here")
                        for y in range(len(self.predictions)):
                            if self.predictions[y][0] == key:
                                # print("here2")

                                if self.gene_annotation[x][3] == self.predictions[y][3]:
                                    start = start+1
                                if self.gene_annotation[x][4] == self.predictions[y][4]:
                                    end = end+1
                                if self.gene_annotation[x][4] == self.predictions[y][4] and self.gene_annotation[x][3] == self.predictions[y][3]:
                                    both = both+1

        print("Matches for : Both {}, Start {}, End {}".format(both, start, end))


# Hidden Markov Model class that uses Viterbi Algorithm and gives you gene predictions given a genome
# (This class uses data from Genome_Data)

class HiddenMarkovModel(object):

    # Class Constructor, initialize states in your HMM

    def __init__(self, states):
        self.states = states
        self.transition_probabilities = dict()
        self.emission_probabilities = dict()

        for state in states:
            probabilities = dict()
            for x in states:
                probabilities.update({x: 0})
            self.transition_probabilities.update({state: probabilities})
            self.emission_probabilities.update({state: dict()})

    # Method to enter transition probability between 2 states

    def enter_transitionProbability(self, State_0, State_1, probability):

        self.transition_probabilities[State_0][State_1] = probability

    # Given average length of intergenic region and genic region, fill out all possible transition probabilities

    def fill_transitionProbabilityTable(self, intergenic_length, genic_length):

        self.transition_probabilities["Intergenic"]["Intergenic"] = float(
            float(intergenic_length-1)/(intergenic_length))
        self.transition_probabilities["Intergenic"]["Start"] = float(
            1.0/(intergenic_length))
        self.transition_probabilities["Start"]["Middle"] = 1.0
        self.transition_probabilities["Middle"]["Middle"] = float(
            float(genic_length-7)/(genic_length-6))
        self.transition_probabilities["Middle"]["Stop"] = float(
            1.0/(genic_length-6))
        self.transition_probabilities["Stop"]["Intergenic"] = 1.0

    # Given frequency tables for integrenic nucleotide, genic start/middle/stop codon, fill out emmission probability table

    def fill_emissionProbabilityTable(self, intergenic_emission, start_emission, middle_emission, stop_emission):

        self.emission_probabilities["Intergenic"] = intergenic_emission
        self.emission_probabilities["Start"] = start_emission
        self.emission_probabilities["Middle"] = middle_emission
        self.emission_probabilities["Stop"] = stop_emission

    # Method for easier access to columns for a 2D data structure

    def column(self, matrix, i):
        return [row[i] for row in matrix]

    # Assuming you have filled out transition probabilities and emission probabilities,
    # This method runs Viterbi Algorithm on a given sequence to fill out a 2D table using Dynammic Programming
    # You can backtrack and retrieve predicted genes in the sequence given

    def Viterbi(self, sequence, initial):

        self.table = np.zeros([len(self.states), len(sequence)])
        self.direction_table = [
            [0 for x in range(len(sequence))] for y in range(len(self.states))]

        for y in range(3):

            for x in range(len(self.column(self.direction_table, y))):
                self.direction_table[x][y] = x
            if y == 0:

                for x in range(len(self.column(self.table, y))):
                    if initial[x]*self.emission_probabilities[self.states[y]][sequence[y]] == 0.0:
                        value = -float("Inf")
                        self.table[x][y] = value
                    else:
                        value = math.log10(
                            initial[y]*self.emission_probabilities[self.states[y]][sequence[0]])
                        self.table[y][0] = value
            else:
                for x in range(len(self.column(self.table, y))):
                    if x == 0:
                        V = self.table[x][y-1]
                        T = math.log10(
                            self.transition_probabilities[self.states[x]][self.states[x]])
                        E = math.log10(
                            self.emission_probabilities[self.states[x]][sequence[y]])

                        value = V+T+E

                        self.table[x][y] = value

                    else:
                        value = -float("Inf")
                        self.table[x][y] = value

        for y in range(3, len(sequence)):

            nucleotide = sequence[y]
            codon = sequence[y-2:y+1]

            for x in range(len(self.column(self.table, y))):
                list = []
                if x == 0:

                    for z in range(len(self.column(self.table, y))):

                        V = (self.table[z][y-1])

                        if self.transition_probabilities[self.states[z]][self.states[x]] == 0.0:
                            T = -float("Inf")
                        else:
                            T = math.log10(
                                self.transition_probabilities[self.states[z]][self.states[x]])

                        if self.emission_probabilities[self.states[x]][nucleotide] == 0.0:
                            E = -float("Inf")
                        else:
                            E = math.log10(
                                self.emission_probabilities[self.states[x]][nucleotide])

                        value = V+T+E

                        list.append(value)

                        max_index = list.index(max(list))

                    self.table[x][y] = max(list)
                    self.direction_table[x][y] = max_index

                else:

                    for z in range(len(self.column(self.table, y))):

                        V = (self.table[z][y-3])

                        if self.transition_probabilities[self.states[z]][self.states[x]] == 0.0:
                            T = -float("Inf")
                        else:
                            T = math.log10(
                                self.transition_probabilities[self.states[z]][self.states[x]])

                        if codon in self.emission_probabilities[self.states[x]]:

                            E = math.log10(
                                self.emission_probabilities[self.states[x]][codon])

                        else:
                            E = -float("Inf")
                    
                        value = V+T+E

                        list.append(value)
                        max_index = list.index(max(list))

                    self.table[x][y] = max(list)
                    self.direction_table[x][y] = max_index

    # This method assumes you have called Viterbi, and filled out dynammic programming table
    # You will backtrack and retrieve the sequence of states in the genome, essentially solcing your hissen markov model for a given state
    # This Method also saves the start and stop position of predicted genes

    def BackTrack(self):

        self.state_sequence = []
        self.start_index = []
        self.stop_index = []
        self.genic_count = 0

        index = (len(self.direction_table[0]))-1

        numcol = (len(self.table[0]))
        numrow = (len(self.table))

        backpointer = self.column(self.table, index).index(
            max(self.column(self.table, index)))

        self.state_sequence = [self.states[backpointer]]+self.state_sequence

        while index > -1:

            state = self.column(self.direction_table, index)[backpointer]

            if self.states[state] == 'Start':

                self.start_index.append(index)

            if self.states[state] == 'Stop':
                self.stop_index.append(index)

            self.state_sequence = [self.states[state]]+self.state_sequence

            if backpointer == 0:
                index = index-1
            else:
                index = index-3

            backpointer = state

        count1 = 0
        count2 = 0

        for x in self.state_sequence:
            if x == "Start":
                count1 = count1+1
            if x == "Stop":
                count2 = count2+1
                self.genic_count = self.genic_count+1

        self.start_index.reverse()
        self.stop_index.reverse()

    # Method that writes predicted genes onto a file in gff3 format, given the sequence and the gene name
    # This Method calls Viterbi and Backtrack

    def writeToFile(self, filename, name, sequence, initial):

        if os.path.exists(filename):
            f = open(filename, "a")
            self.Viterbi(sequence, initial)
            self.BackTrack()
            for x in range(len(self.stop_index)):

                f.write("{} ena CDS {} {} . + 0 . \n".format(name,
                                                             self.start_index[x], self.stop_index[x]))

        else:
            f = open(filename, "w+")
            self.Viterbi(sequence, initial)
            self.BackTrack()
            for x in range(len(self.stop_index)):

                f.write("{} ena CDS {} {} . + 0 . \n".format(name,
                                                             self.start_index[x], self.stop_index[x]))

        f.close()

# Main method to answer questions on Assignment # 3 for COMP 462-Computational Biology Methods

def main():

    # QUESTION 1a): Getting Average Length for Intergenic & Genic Regions, Nucleotide Frequency Table, & Codon Frequency Table
    #             given fasta file(genome) and gff3 file(gene annotation)

    data = Genome_Data("Vibrio_cholerae.GFC_11.dna.toplevel.fa",
                       "Vibrio_cholerae.GFC_11.37.gff3")

    data.Integernetic_and_Genetic_Region()
    data.printAverageLengths()
    print("\n")

    data.updateFrequencyTables(data.intergenic_region, data.genic_region)
    print("Nucleotide Frequency Table for Intergenic Region\n")
    data.printNucleotide_FrequencyTable()
    print("\n")

    print("Codon Frequency Table for Genic Region\n")
    data.printCodon_FrequencyTable()
    print("\n")

    print("Start Codon Frequency Table for Genic Region\n")
    data.printStartCodonFrequency()
    print("\n")

    print("Middle Codon Frequency Table for Genic Region\n")
    data.printMiddleCodonFrequency()
    print("\n")

    print("Stop Codon Frequency Table for Genic Region\n")
    data.printStopCodonFrequency()
    print("\n")


    # QUESTION 1b): Implement Viterbi Algorithm with a HMM with 4 states as seen in class
    #               write a GFF3 file with one line per predicted gene


    hmm = HiddenMarkovModel(["Intergenic", "Start", "Middle", "Stop"])

    hmm.fill_transitionProbabilityTable(
        data.intergenic_length, data.genic_length)

    hmm.fill_emissionProbabilityTable(
        data.nucleotide_frequency, data.startCodon_frequency, data.middleCodon_frequency, data.stopCodon_frequency)

    
    for keys in data.seq_dict:
        hmm.writeToFile("Predicted_1B.gff3",keys,data.seq_dict[keys],[1,0,0,0])
    

    # QUESTION 1c): Predict genes for Vibrio vulnificus  using an HMM with the parameters obtained for Vibrio cholerae
    #               you will need contruct another Genome_Data object and get your data and run HMM on it
    #               you will also get a gff3 file of predicted genes

    data2 = Genome_Data(
        "Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa", "Vibrio_cholerae.GFC_11.37.gff3")

    for keys in data2.seq_dict:
        hmm.writToFile("Predicted_1C.gff3",keys,data2.seq_dict[keys],[1,0,0,0])
    
    

    # QUESTION 1d): Given the actual gene annotation for Vibrio vulnificus, compare how your predictions do
    #               compared to the real file, compare start/end of the predicted and actual genes
    #               you will need to create another Genome_Data object again to save new gene annotations and use the
    #               compareWithPrediction( method)


    data3 = Genome_Data("Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa",
                        "Vibrio_vulnificus.ASM74310v1.37.gff3")

    data3.compareWithPrediction("Predicted_1c.gff3")

    print("\n")
    print("Gene Annotations: {}, Gene Predictions: {}".format(len(data3.gene_annotation),len(data3.predictions)))
    print("\n")


    # QUESTION 1e): What properties make genes easier to predict than others, and some easier to predict using Viterbi Algorithm
    #               to answer this question we will look at start and stop codons for given fasta file and using gene annotation &
    #               gene prediction, and generate frequency tables to identify factors in gene prediction

    

    data3.Integernetic_and_Genetic_Region()
    data3.updateFrequencyTables(data.intergenic_region, data.genic_region)

    print("Nucleotide Frequency Table for Annotated Gene Intergenic\n")
    data3.printNucleotide_FrequencyTable()
    print("\n")

    print("Start Codon Frequency Table for Annotated Gene\n")
    data3.printStartCodonFrequency()
    print("\n")

    print("Middle Codon Frequency Table for Annotated Gene\n")
    data3.printMiddleCodonFrequency()
    print("\n")

    print("Stop Codon Frequency Table for Annotated Gene\n")
    data3.printStopCodonFrequency()
    print("\n")

    data4=Genome_Data("Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa","Predicted_1c.gff3")
    data4.Integernetic_and_Genetic_Region()
    data4.updateFrequencyTables(data.intergenic_region, data.genic_region)

    print("Nucleotide Frequency Table for Predicted Gene Intergenic\n")
    data4.printNucleotide_FrequencyTable()
    print("\n")

    print("Start Codon Frequency Table for Predicted Gene\n")
    data4.printStartCodonFrequency()
    print("\n")

    print("Middle Codon Frequency Table for Predicted Gene\n")
    data4.printMiddleCodonFrequency()
    print("\n")

    print("Stop Codon Frequency Table for Predicted Gene\n")
    data4.printStopCodonFrequency()
    print("\n")




if __name__ == '__main__':
    main()
