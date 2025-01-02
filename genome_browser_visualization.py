import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse
import matplotlib.patheffects as patheffects
import matplotlib.image as mplimg


parser = argparse.ArgumentParser()

parser.add_argument('--outputFile','-o',default='genome_visualization.png',type=str,action='store',help='output file goes here')
parser.add_argument('--inputFile','-i',default='/path/toyour/data',type=str,action='store',help='input file goes here')
parser.add_argument('--outTextFile','-t',default='Lecture4.data',type=str,action='store',help='output file goes here')


args = parser.parse_args()

outFile=args.outputFile
inFile=args.inputFile
outText=args.outTextFile

print(outFile,inFile)

iBlue=(88/255,85/255,120/255)

figureWidth=5
figureHeight=5

plt.style.use('BME163')

plt.figure(figsize=(figureWidth,figureHeight))

panelWidth = 4
panelHeight = 1.5

panel1 = plt.axes([0.1/figureWidth,3.3/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])
panel2 = plt.axes([0.1/figureWidth,1.7/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])
panel3 = plt.axes([0.1/figureWidth,0.1/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])

panel1.set_xticks([])
panel1.set_yticks([])
panel2.set_xticks([])
panel2.set_yticks([])
panel3.set_xticks([])
panel3.set_yticks([])

#####PANEL2#####PSLFILE#####
panel2.set_xlim(45232000, 45241000)
panel2.set_ylim(-1, 70)

psl_vals = []
cov_dict={}

#parse through psl file
with open('BME163_Input_Data_6.psl', 'r') as psl_file:
    for i, line in enumerate(psl_file):
        line = line.split('\t')  

        #relevant data
        chromosome = line[13]
        start = int(line[15])
        end = int(line[16])
        blockstarts = np.array(line[20].split(',')[:-1], dtype=int)
        blockwidths = np.array(line[18].split(',')[:-1], dtype=int)

        ##also add ranges for start and instead instead of just chr7
        if chromosome == 'chr7' and (start >= 45232000 and start <= 45241000) or (end >= 45232000 and end <= 45241000):
            psl_vals.append([start, end, blockstarts, blockwidths, False]) #BOOLEAN false for stacking later

#sort psl values
sorted_psl = sorted(psl_vals,key = lambda x:x[1])

cov_dict={} #for the histogram later

#iterate through sorted psl values list
for y in range(0, len(sorted_psl), 1):
    separate = 0
    ymax = y
    for val in sorted_psl:
        left = val[0]
        right = val[1]
        if not val[4]:
            if left > separate:
                height = 0.05 #line connecting rectangles
                rectangle = mplpatches.Rectangle([left, y-(height/2)], 
                                                 right-left, 
                                                 height, 
                                                 facecolor = 'black',
                                                 linewidth=0)
                panel2.add_patch(rectangle)

                #blocks
                height=0.50
                for index in range(0, len(val[2]), 1):
                    blockstart = val[2][index]
                    blockwidth = val[3][index]
                    blockend = blockstart + blockwidth
                    rectangle = mplpatches.Rectangle([blockstart, 
                                                      y-(height/2)], 
                                                      blockwidth, 
                                                      height, 
                                                      facecolor=iBlue, 
                                                      edgecolor='black', 
                                                      linewidth=0.1)
                    panel2.add_patch(rectangle)
                    separate = right
                    val[4]=True

                    #add to dictionary for the histogram
                    for nuc in range(blockstart, blockend, 1):
                        if nuc in cov_dict:
                            cov_dict[nuc]+=1
                        else:
                            cov_dict[nuc]=1



#####PANEL3#####PSLFILEHISTOGRAM#####
panel3.set_xlim(45232000, 45241000)
panel3.set_ylim(0,61)

#same as code for assignment 2, iterate differently
for nucs, counts in cov_dict.items():
    bot=0
    height=counts
    width=1
    rectangle=mplpatches.Rectangle([nucs,bot],
                                   width, 
                                   height, 
                                   facecolor=iBlue, 
                                   edgecolor=iBlue, 
                                   linewidth=0)
    panel3.add_patch(rectangle)


                
#####PANEL1#####GTFFILE#####
#basically just modified the psl code A LOT but same thinking
panel1.set_xlim(45232000, 45241000)
panel1.set_ylim(-1,10)

gtf_vals = []
transcript = []

#parse through gtf file
with open('gencode.vM12.annotation.gtf', 'r') as gtf_file:
    for _ in range(5):  
        next(gtf_file)

    
    for line in gtf_file:
        line = line.split('\t')
        chromosome = line[0]
        feature = line[2]

        #relevant data, could've put before 
        if feature in ['exon', 'CDS'] and chromosome == 'chr7':
            start = int(line[3])
            end = int(line[4])
            blockstart = start
            blockwidth = end - start
            metaData = line[8].split('transcript_id "')[1].split('"')[0] #lecture 23 code

            #2d list
            gtf_vals.append((chromosome, start, end, blockstart, blockwidth, feature, metaData))

            found = False #Boolean for stacking
            for sub in transcript:
                if sub[0] == metaData:
                    sub[2] = (min(sub[2][0], start), max(sub[2][1], end))
                    #Exon storage
                    if feature == 'exon':
                        sub[3][0].append(blockstart)
                        sub[3][1].append(blockwidth)
                        sub[3][2].append(0.25)
                    #CDS storage
                    elif feature == 'CDS':
                        sub[4][0].append(blockstart)
                        sub[4][1].append(blockwidth)
                        sub[4][2].append(0.5)
                    sub[1] = False
                    found = True
                    break
            
            #could optimize data structurem but lists inside of lists and add to trasncript list
            if not found:
                if feature == 'exon':
                    sub = [metaData, False, (start, end), [[blockstart], [blockwidth], [0.25]], [[], [], []]]
                elif feature == 'CDS':
                    sub = [metaData, False, (start, end), [[], [], []], [[blockstart], [blockwidth], [0.5]]]
                transcript.append(sub)


#iterate through transcript list, dont need to sort
for y in range(0, len(transcript), 1):
    separate = 0
    ymax = y
    
    for val in transcript:
        left = val[2][0]
        right = val[2][1]
        if not val[1]:
            if left > separate:
                
                #line connecting
                height = 0.05
                rectangle = mplpatches.Rectangle([left, y-(height/2)], 
                                                 right-left, 
                                                 height, 
                                                 facecolor = 'Grey',
                                                 edgecolor='black', 
                                                 linewidth=0.25)
                panel1.add_patch(rectangle)

                # Exon plots
                for j in range(0, len(val[3][0]), 1):
                    start = val[3][0][j]
                    width = val[3][1][j]
                    end = start + width
                    height = val[3][2][j]
                    rectangle = mplpatches.Rectangle([start, y-(height/2)],
                                                      width, 
                                                      height, 
                                                      facecolor='Grey', 
                                                      edgecolor='black', 
                                                      linewidth=0.25)
                    panel1.add_patch(rectangle)

                # CDS plots
                for k in range(0,len(val[4][0]),1):
                    start = val[4][0][k]
                    width = val[4][1][k]
                    end = start + width
                    height = val[4][2][k]
                    rectangle = mplpatches.Rectangle([start, y-(height/2)], 
                                                     width, 
                                                     height, 
                                                     facecolor='Grey', 
                                                     edgecolor='black', 
                                                     linewidth=0.25)
                    panel1.add_patch(rectangle)

                separate = right
                val[1]=True



plt.savefig(outFile,dpi=2400)
