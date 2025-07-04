# Generating a textfile for uploading bw files to the UCSC browser

string1 = """track type=bigWig name="track_title" bigDataUrl=https://storage.googleapis.com/siwase/rai1ko_hesc/sample visibility=full color=147,112,219 autoScale=on alwaysZero=off gridDefault=off maxHeightPixels=60 graphType=bar viewLimits=0:5 windowingFunction=mean smoothingWindow=2\n\n"""

x=[
'Rai1-het_d0_rep1.bw',
'Rai1-het_d0_rep2.bw',
'Rai1-het_d0_rep3.bw',
'Rai1-het_d0_rep4.bw',
'Rai1-het_d20_rep1.bw',
'Rai1-het_d20_rep2.bw',
'Rai1-het_d20_rep3.bw',
'Rai1-het_d20_rep4.bw',
'Rai1-het_d34_rep1.bw',
'Rai1-het_d34_rep2.bw',
'Rai1-het_d34_rep3.bw',
'Rai1-het_d34_rep4.bw',
'Rai1-het_d60_rep1.bw',
'Rai1-het_d60_rep2.bw',
'Rai1-het_d60_rep3.bw',
'Rai1-het_d60_rep4.bw',
'Rai1-ko_d0_rep1.bw',
'Rai1-ko_d0_rep2.bw',
'Rai1-ko_d0_rep3.bw',
'Rai1-ko_d0_rep4.bw',
'Rai1-ko_d20_rep1.bw',
'Rai1-ko_d20_rep2.bw',
'Rai1-ko_d20_rep3.bw',
'Rai1-ko_d20_rep4.bw',
'Rai1-ko_d34_rep1.bw',
'Rai1-ko_d34_rep2.bw',
'Rai1-ko_d34_rep3.bw',
'Rai1-ko_d34_rep4.bw',
'Rai1-ko_d60_rep1.bw',
'Rai1-ko_d60_rep2.bw',
'Rai1-ko_d60_rep3.bw',
'Rai1-ko_d60_rep4.bw',
'Rai1-wt_d0_rep1.bw',
'Rai1-wt_d0_rep2.bw',
'Rai1-wt_d0_rep3.bw',
'Rai1-wt_d0_rep4.bw',
'Rai1-wt_d20_rep1.bw',
'Rai1-wt_d20_rep2.bw',
'Rai1-wt_d20_rep3.bw',
'Rai1-wt_d20_rep4.bw',
'Rai1-wt_d34_rep1.bw',
'Rai1-wt_d34_rep2.bw',
'Rai1-wt_d34_rep3.bw',
'Rai1-wt_d34_rep4.bw',
'Rai1-wt_d60_rep1.bw',
'Rai1-wt_d60_rep2.bw',
'Rai1-wt_d60_rep3.bw',
'Rai1-wt_d60_rep4.bw'
]

with open("upload.txt", "w") as text_file:
 for i in x:
  string2 = string1.replace("track_title", "%s" % (i))
  string3 = string2.replace("sample", "%s" % (i))
  text_file.write(string3)
