import numpy as np
import math
import sys

## Needleman Wunsch Dynammic Programming Algorithm for Sequence Alignment under Affine Gap Penalty
## Code used to Answer Question 3C

class NW_AffineGap(object):
    def __init__(self,match,mismatch,gap_penalty1,gap_penalty2):
        self.match=match
        self.mismatch=mismatch
        self.gap_penalty1=gap_penalty1
        self.gap_penalty2=gap_penalty2
    
    def initializeTables(self,sequence1, sequence2):
        self.track=[[0 for x in range(len(sequence1)+1)] for y in range(len(sequence2)+1)]
        self.table=np.zeros([len(sequence2)+1,len(sequence1)+1])
        self.tableIx=np.zeros([len(sequence2)+1,len(sequence1)+1])
        self.tableIy=np.zeros([len(sequence2)+1,len(sequence1)+1])
        self.sequence1=sequence1
        self.sequence2=sequence2

        for x in range(len(sequence2)):
    
            self.tableIx[x][0]=self.gap_penalty1+(x*self.gap_penalty2)
            self.table[x+1][0]=-float("inf")
            self.tableIy[x+1][0] =-float("inf")

            self.track[x+1][0]="up"

        self.tableIx[len(sequence2)][0]=self.gap_penalty1+(len(sequence2)*self.gap_penalty2)

        for y in range(len(sequence1)):
            
            self.tableIy[0][y]=self.gap_penalty1+(y*self.gap_penalty2)
            self.table[0][y+1]=-float("inf")
            self.tableIx[0][y+1]=-float("inf")
            self.track[0][y+1]="left"

        self.tableIy[0][len(sequence1)]=self.gap_penalty1+(len(sequence1)*self.gap_penalty2)

        
    def match_or_mismatch(self,basepair1,basepair2):
        if basepair1==basepair2:
            return self.match
        else:
            return self.mismatch   


    
    def fillTable(self):

        for x in range(1,len(self.sequence2)+1):
            for y in range(1,len(self.sequence1)+1):         

                

                alignments=[self.table[x-1][y-1]+self.match_or_mismatch(self.sequence2[x-1],self.sequence1[y-1]),
                self.tableIx[x-1][y-1]+self.match_or_mismatch(self.sequence2[x-1],self.sequence1[y-1]),
                self.tableIy[x-1][y-1]+self.match_or_mismatch(self.sequence2[x-1],self.sequence1[y-1])]

                self.table[x][y]=max(alignments)

                self.tableIx[x][y]=max([self.table[x-1][y]+self.gap_penalty1,self.tableIx[x-1][y]+self.gap_penalty2])
                self.tableIy[x][y]=max([self.table[x][y-1]+self.gap_penalty1,self.tableIy[x][y-1]+self.gap_penalty2])

                

                if(alignments.index(max(alignments))==0):
                    self.track[x][y]="diagonal"
                elif(alignments.index(max(alignments))==1):
                    self.track[x][y]="up"
                elif(alignments.index(max(alignments))==2):
                    self.track[x][y]="left"
               
            
        self.track=np.array(self.track)
            
                    

    def backTrack(self):
        alignment1=""
        alignment2=""

        
        x=len(self.track)-1
        y=len(self.track[0])-1

        while (1):
            
            if(self.track[x][y]=="diagonal"):
                alignment1=self.sequence1[y-1]+alignment1
                alignment2=self.sequence2[x-1]+alignment2

                x=x-1
                y=y-1
                if(x==0 | y==0):
                    break

            elif(self.track[x][y]=="up"):
                if(self.track[x-1][y]=="up"):
                    self.track[x-1][y]="diagonal"

                alignment1="-"+alignment1
                alignment2=self.sequence2[x-1]+alignment2

                x=x-1
                y=y
                if(x==0 | y==0):
                    break

            elif(self.track[x][y]=="left"):
                if(self.track[x][y-1]=="left"):
                    self.track[x][y-1]="diagonal"

                alignment1=self.sequence1[y-1]+alignment1
                alignment2="-"+alignment2

                x=x
                y=y-1

                if(x==0 | y==0):
                    break
                

        print(alignment1)
        self.alignment1=alignment1
        print(alignment2)
        self.alignment2=alignment2


    def scoreAlignments(self):
        score=0
        for x in range(len(self.alignment1)):
            if(self.alignment1[x]==self.alignment2[x]):
                score=score+self.match
            else:
                score=score+self.mismatch
        return score

       

        
def main():

    sequenceAlignment=NW_AffineGap(1,-1,-1,-float("inf"))
    f = open("hw1_long.fa","r")
    sequences=[]
    for word in f:
        if(word[0]!='>'):
            sequences.append(word)
    
    sequenceAlignment.initializeTables(sequences[0],sequences[2])
    sequenceAlignment.fillTable()
    sequenceAlignment.backTrack()
 
    print(sequenceAlignment.scoreAlignments())

   


if __name__=='__main__':
    main()
