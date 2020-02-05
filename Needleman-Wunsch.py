import numpy as np

## Needleman Wunsch Dynammic Programming Algorithm for Sequence Alignment under Linear Gap Penalty

class Needleman_Wunsch(object):
    def __init__(self,match,mismatch,gap_penalty):
        self.match=match
        self.mismatch=mismatch
        self.gap_penalty=gap_penalty
    
    def initializeTables(self,sequence1, sequence2):
        self.track=[[0 for x in range(len(sequence1)+1)] for y in range(len(sequence2)+1)]
        self.table=np.zeros([len(sequence2)+1,len(sequence1)+1])
        self.sequence1=sequence1
        self.sequence2=sequence2

        for x in range(len(sequence2)):
            self.table[x+1][0]=-(x+1)
            self.track[x+1][0]="up"

        for y in range(len(sequence1)):
            self.table[0][y+1]=-(y+1)
            self.track[0][y+1]="left"

    def match_or_mismatch(self,basepair1,basepair2):
        if basepair1==basepair2:
            return self.match
        else:
            return self.mismatch    
    
    def fillTable(self):

        for x in range(1,len(self.sequence2)+1):
            for y in range(1,len(self.sequence1)+1):

                alignments=[self.table[x-1][y-1]+self.match_or_mismatch(self.sequence2[x-1],self.sequence1[y-1]),
                self.table[x-1][y]+self.gap_penalty,
                self.table[x][y-1]+self.gap_penalty]

                self.table[x][y]=max(alignments)

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
                alignment1="-"+alignment1
                alignment2=self.sequence2[x-1]+alignment2

                x=x-1
                y=y
                if(x==0 | y==0):
                    break

            elif(self.track[x][y]=="left"):
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

    sequenceAlignment=Needleman_Wunsch(1,-1,-1)
    #f = open("hw1_medium.fa","r")
    #sequences=[]
    #for word in f:
        #if(word[0]!='>'):
            #sequences.append(word)
    
    sequenceAlignment.initializeTables("HAPE","APPLE")
    sequenceAlignment.fillTable()
    print(sequenceAlignment.table)
    print(sequenceAlignment.track)
    sequenceAlignment.backTrack()
    print(sequenceAlignment.scoreAlignments())

    
   

if __name__=='__main__':
    main()
