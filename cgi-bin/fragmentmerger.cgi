#!/usr/bin/python -u

# Output buffering disabled by "-u" above

import datetime

starttime = datetime.datetime.now()

import cgi, sys, socket, re, os, time

import climb

# import cgitb
# cgitb.enable()
# cgitb.enable(display=1, logdir='/var/www/tmp')

# sys.stderr = sys.stdout

form = cgi.FieldStorage()

version = "Fragment Merger-DEV"

def closePage():
	print '<a href="/fmt/index.html">Submit another</a>'
	print '<hr>'
	print '<tt><p align=right>Version 1.0 (18 January 2011)<br>'
	endtime = datetime.datetime.now()
	print 'Request from %s<br>' % (cgi.os.environ['REMOTE_ADDR'])
	print 'Served by %s at %s<br>' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0])
	print 'Run completed %s<br>' % (endtime.strftime("%Y-%m-%d %H:%M"))
 	print '%03.6f seconds</tt></p>' % ((endtime - starttime).microseconds / float(1000000))
	print '</body></html>'

	f = open('/var/www/tmp/fm.log', 'a')
	f.write( 'Fragment Merger-DEV\t' )
	f.write( 'Request from %s ' % (cgi.os.environ['REMOTE_ADDR']) )
	f.write( 'Random Token %s\t' % (randomToken) )
	f.write( 'Served by %s at %s\t' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0]) )
	f.write( 'Run completed %s\t' % (endtime.strftime("%Y-%m-%d %H:%M")) )
	f.write( '%03.6f seconds\t' % ((endtime - starttime).microseconds / float(1000000)) )
	f.write( '\n')
	f.close()

def badFile(errorText, explainText):
	print "<h1>%s</h1>" % (version)
	print "<h2>%s</h2>" % (errorText)
	print explainText + '<br>'
	closePage()
	sys.exit()

def htmlPrint(pp, lines=1):
	print pp
	for i in range(lines):
		print "<br>"

print 'Content-Type: text/html'
print

print '<!DOCTYPE html'
print '	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"'
print '	 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">'
print '<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">'
print '<head><title>HVDR %s</title>' % (version)
print '<link rel="stylesheet" type="text/css" href="/res/hvdr.css"/>'
print '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/>'

print '<script type="text/javascript">'
print '// drawPercentBar()'
print '// Written by Matthew Harvey (matt AT unlikelywords DOT com)'
print '// (http://www.unlikelywords.com/html-morsels/)'
print '// Distributed under the Creative Commons'
print '// "Attribution-NonCommercial-ShareAlike 2.0" License'
print '// (http://creativecommons.org/licenses/by-nc-sa/2.0/)'
print 'function drawPercentBar(width, percent, color, background)'
print '{'
print '  var pixels = width * (percent / 100);'
print '  if (!background) { background = "none"; }'
print '  document.write("<div style=\\"position: relative; line-height: 1em; background-color: " + background + "; border: 1px solid black; width: " + width + "px\\">");'
print '  document.write("<div style=\\"height: 1.5em; width: " + pixels + "px; background-color: " + color + ";\\"></div>");'
print '  document.write("<div style=\\"position: absolute; text-align: center; padding-top: .25em; width: " + width + "px; top: 0; left: 0\\">" + percent + "%</div>");'
print '  document.write("</div>");'
print '}'
print '</script>'

# --------- jQuery ----------
print '<script src="/jquery-1.5.2.min.js"></script>'
print '<script>'
print '$(document).ready(function()'
print '{'
print '    $("div#paragraph1").hide();'
print '    $("a#paragraph1Anchor").click(function()'
print '    {'
print '        $("div#paragraph1").toggle(\'fast\');'
print '        return false;'
print '    });'
print '    $("div#paragraph2").hide();'
print '    $("a#paragraph2Anchor").click(function()'
print '    {'
print '        $("div#paragraph2").toggle(\'fast\');'
print '        return false;'
print '    });'
print '    $("div#paragraph3").hide();'
print '    $("a#paragraph3Anchor").click(function()'
print '    {'
print '        $("div#paragraph3").toggle(\'fast\');'
print '        return false;'
print '    });'
print '    $("div#paragraph4").hide();'
print '    $("a#paragraph4Anchor").click(function()'
print '    {'
print '        $("div#paragraph4").toggle(\'fast\');'
print '        return false;'
print '    });'
print '});'

print '</script>'
# ---------------------------

print '</head><body>'

# os.system('beep -f 4096 -r 2 -l 32 -d 64')	# make a noise; was -r 16

randomToken = climb.randomStamp()
climb.TEMPFOLDER = '/var/www/tmp/'

shortAmplicon = int(form["shortamplicon"].value)
trimPercentage = float(form["trimmedlengththreshold"].value)

mergedID = form["mergedID"].value
referenceSequenceFile = form["referencesequencefile"].filename

fragmentFilenames = []
for i in range(1,climb.MAX_FRAGMENTS+1):
	fragmentFilenames.append(form["fragmentfile%i" % i].filename)

lastFragment = 0	# index position -- zero-indexed
for i in range(climb.MAX_FRAGMENTS):
	if len(fragmentFilenames[i]) == 0:
		lastFragment = i-1		# Previous fragment considered the last; negative i handled below
		break

if (len(fragmentFilenames[0]) == 0 and len(fragmentFilenames[1]) == 0) or (lastFragment < 1):
	badFile('Specify at least the first two fragment files!', "")

fragmentFile = []
for i in range(1,climb.MAX_FRAGMENTS+1):
	fragmentFile.append(form["fragmentfile%i" % i].file)

htmlPrint('<h1>HVDR Fragment Merger Tool</h1>', lines=0)
htmlPrint('<h2><img src="/res/75px-Symbol_merge_review.png" style="vertical-align:middle"> Fragment Merges Completed</h2>')

# chromatogram file is stored locally so that the filename can be passed to the Sequence.load() method
# all files stored in a folder with a randomly-generated name
# ensure chromatograms and alignment temp files are written to session folder in climb.contig() method
climb.TEMPFOLDER += 'Fragment-Merger-' + randomToken + '/'
os.makedirs(climb.TEMPFOLDER)

tempF = []
for i in range(lastFragment + 1):
	tempFF = climb.TEMPFOLDER + fragmentFilenames[i]
	tempF.append(tempFF)	# required later by the merge method
	temp = open(tempFF, 'wb')
	temp.write(fragmentFile[i].read())
	temp.close()
	fragmentFile[i].close()

if len(referenceSequenceFile) > 0:
	refSeqFilename = climb.TEMPFOLDER + referenceSequenceFile
	temp = open(refSeqFilename, 'wb')
	temp.write(form["referencesequencefile"].file.read())
	temp.close()

	R = climb.Sequence()
	R.load(climb.TEMPFOLDER + referenceSequenceFile)	# will terminate if not a FASTA file; crude but effective
	R.unload()						# if FASTA file, then unload

slideM = form['slideMotif'].value.upper()
slideP = form['slidePosition'].value
if slideP.isdigit():
	slideP = int(slideP)
else:
	slideP = 0

# Populate lists with reverse (true/false) and compliment (true/false)
Rev = [form.has_key('rev%i' % i) for i in range(1,climb.MAX_FRAGMENTS+1)]
Comp = [form.has_key('comp%i' % i) for i in range(1,climb.MAX_FRAGMENTS+1)]

# List containing TRUE if fragment is ABIF, else FALSE
fragABIF = [form["type%i" % i].value == 'ABIF' for i in range(1,climb.MAX_FRAGMENTS+1)]

# Populate list if trim window and trim threshold values
TW = [int(form["tw%i" % i].value) for i in range(1,climb.MAX_FRAGMENTS+1)]
TT = [int(form["tt%i" % i].value) for i in range(1,climb.MAX_FRAGMENTS+1)]

MergeContainer = climb.fragmentmerger(tempF[:lastFragment+1], FRAGABIF=fragABIF, R=Rev, C=Comp, X=TW, Y=TT, mergeID=mergedID, slideMotif=slideM, slidePosition=slideP)
# if mergeID his given a value, it will be used, otherwise, the value from the file will be used
Merge = MergeContainer[0]	# [0] is merged sequence object, [1], [2] and [3] are fragment objects

# ID of first sequence modified in method to contain all three IDs; megamerger takes ID from first fragment
# Using fragmentFile1filename instead
# optional mergedID used if specified
if len(mergedID) > 0:
	mergedID = mergedID.strip()		# remove laeding and trailing whitespace
	mergedID = mergedID.replace(' ', '_')	# replace internal whitespace with underscore
	mergedID = mergedID.replace('-', '_')	# replace internal 'gaps' with underscore
	Merge.seq[0]['id'] = mergedID		# changed in "fragmentermerger" method which writes the data out

htmlPrint('Hold the mouse over the icons below for an explanation of each. More details are provided below the table. The download and output links provided below the table will expire after one hour.', lines=2)

print('<table border="1">')

for i in range(1,lastFragment+2):
	if i % 2 == 0:
		rowCol = 'plum'
	else:
		rowCol = 'lightsteelblue'

	print('<tr bgcolor="%s"><td></td><td><b>Fragment %i</b></td><td>%s</td></tr>' % (rowCol, i, fragmentFilenames[i-1]))

	print('<tr bgcolor="%s"><td></td><td>Fragment Type</td><td>%s</td>' % (rowCol, form["type%i" % i].value))

	if Rev[i-1] == True:
		revOut = '<b>Yes</b>'
	else:
		revOut = 'No'

	if Comp[i-1] == True:
		compOut = '<b>Yes</b>'
	else:
		compOut = 'No'

	print('<tr bgcolor="%s"><td></td><td>Reversed</td><td>%s</td></tr>' % (rowCol, revOut))
	print('<tr bgcolor="%s"><td></td><td>Complemented</td><td>%s</td></tr>' % (rowCol, compOut))

	print('<tr bgcolor="%s"><td></td><td>Original Length (Bases)</td><td>%i</td></tr>' % (rowCol, MergeContainer[i].originalLength))

	# if fragABIF[i-1]:	# only show trim details if file is ABIF
	if True:
		trimmedLengthPercentage = len(MergeContainer[i].seq[0]['seq'])/float(MergeContainer[i].originalLength)*100

		trimmedLength = len(MergeContainer[i].seq[0]['seq'])

		if not fragABIF[i-1]:		# if FASTA file
			trimIcon = '<img src="/res/info.png" alt="info" title="FASTA input file (no trimming)" />'
			barCol = "dodgerblue"
		elif trimmedLengthPercentage < trimPercentage:
			trimIcon = '<img src="/res/warning.png" alt="warning" title="Trimmed length is < %03.2f%% of original length" />' % (trimPercentage)
			barCol = 'red'
		elif trimmedLength >= shortAmplicon:
			trimIcon = '<img src="/res/success.png" alt="success" title="Trimmed length is >= %03.2f%% of original length" />' % (trimPercentage)
			barCol = 'green'
		else:
			trimIcon = '<img src="/res/query.png" alt="query" title="Original length < %i bases" />' % (shortAmplicon)
			barCol = 'yellow'

		print('<tr bgcolor="%s"><td>%s</td><td>Trimmed Length (Bases)</td><td>%i (%3.2f%% of original)<script type="text/javascript">drawPercentBar(200, %i, "%s");</script></td></tr>' % (rowCol, trimIcon, len(MergeContainer[i].seq[0]['seq']), trimmedLengthPercentage, round(trimmedLengthPercentage), barCol))
		print('<tr bgcolor="%s"><td></td><td>Bases trimmed from start of sequence</td><td>%i</td></tr>' % (rowCol, MergeContainer[i].trimLeft))
		print('<tr bgcolor="%s"><td></td><td>Bases trimmed from end of sequence</td><td>%i</td></tr>' % (rowCol, MergeContainer[i].trimRight))

print('</table><br>')
print '<style type="text/css" rel="stylesheet">'
print '<!--'
print '.nowrap {font-family: monospace; font-size: 12px; font-weight: normal; background: #000000; color: #FFFFFF; white-space: nowrap; border:1px solid black; overflow-y:hidden; overflow-x:scroll;}'
print '.output {background: #C1FDC1; margin: 16px; padding: 8px; border: solid; }'
print 'td.Correct {background: green; }'
print 'td.Check {background: red; }'
print 'table.fixed {table-layout:fixed; }'
# http://stackoverflow.com/questions/4185814/fixed-table-cell-width
print '-->'
print '</style>'

############################
# Post-merginging actions: split sequence in the middle, merge fragements, slide
# awk '{ mylen = length($0); first = substr($0, 1, mylen/2); second = substr($0, mylen/2+1); print first, second }' dummy | awk '{ print $1,$2 }'
#
# mergedFile1 = file(climb.TEMPFOLDER + 'mergedFilePart1.fasta', 'w')
# mergedFile2 = file(climb.TEMPFOLDER + 'mergedFilePart2.fasta', 'w')
# mergedfile1.writeline(Merge.seq[0]['id'] + '\n')
# mergedfile2.writeline(Merge.seq[0]['id'] + '\n')
#
# stopped working on this -- the 6th fragment is reversed so there is no 'final' overlap, so this is not necessary
# could be added in future
############################

print '<h3>Merged Sequence</h3>'
print '<a href="#" id="paragraph1Anchor">Click to show/hide merged sequence</a>'
print '<div id="paragraph1" class="output">'
print '<p class="nowrap">'
htmlPrint('>' + Merge.seq[0]['id'])
htmlPrint(Merge.seq[0]['seq'])
print '</p><br>'

htmlPrint('Length of Merged Sequence: %i' % len(Merge.seq[0]['seq']), 2)
if Merge.mergeSlide:
	htmlPrint("Slide successful (Motif %s to position %i)." % (slideM, slideP), 2)
elif len(slideM) > 0:
	htmlPrint("Slide not successful (Motif %s to position %i)." % (slideM, slideP), 2)
else:
	htmlPrint("No slide specified.", 2)
saveFilename = '/tmp/Fragment-Merger-' + randomToken + '/' + os.path.basename(Merge.fileName)
htmlPrint('<a href="%s"><img src="/res/download.png" style="border-style: none"></a> Right-click icon and select "Save (Link) As" to download Merge in FASTA format.' % (saveFilename))

print '</div><br>'

print '<h3>Detailed Merge Output</h3>'
print '<a href="#" id="paragraph2Anchor">Click to show/hide detailed merge output</a>'
print '<div id="paragraph2" class="output">'

# climb.TEMPFOLDER is /var/www/tmp, but we only want /tmp for the actual link...
actualTEMPFOLDER = '/tmp/Fragment-Merger-' + randomToken + '/'
print '<table border="1" class="fixed" width="100%"><col width="64%" col width="12%" col width="12%" col width="12%"><tr><td><em>Merge Alignment</em></td><td><em>Orientation</em></td><td><em>Score</em></td><td><em>Detailed output</em></td></tr>'

for i in range(lastFragment):
	tempFilename = '/var/www' + actualTEMPFOLDER
	tempFilename += 'outFile%02i.txt' % (i+1)
	theFile = open(tempFilename, 'r')
	score = 'N/A'
	while True:
		l = theFile.readline()
		if l[:9] == '# Score: ':			# or re.match("^# Score: ", l)
			score = l[9:].strip()			# treat as string
		if (l[0] != '#') and (len(l) > 1):	# alignment data is first non-blank line not starting with #
			firstLine = l
			firstLineOut = firstLine.replace(' ', '&nbsp;')		# first and second line scanned below to determine merge quality
			alignmentLine = theFile.readline().replace(' ', '&nbsp;')
			secondLine = theFile.readline()
			secondLineOut = secondLine.replace(' ', '&nbsp;')
			break

	theFile.close()

	# first 19 characters are the FASTA ID; strip these out (removes the possibility of a - in the FASTA ID)
	# what remains is then a digit (position/co-ordinate), a space, and the sequence data
	# strip out the digit (1) and the space; therefore, first character should be either a gap or a base

	nonGapFirst  = re.search('[^-]', firstLine[21:])
	if nonGapFirst != None:
		nonGapFirst = nonGapFirst.start()
	else:
		nonGapFirst = -1				# if gap not found
	nonGapSecond = re.search('[^-]', secondLine[21:])
	if nonGapSecond != None:
		nonGapSecond = nonGapSecond.start()
	else:
		nonGapSecond = -1				# if gap not found

	if nonGapFirst < nonGapSecond:
		mergeQuality = 'Correct'
	else:					# equal (even if both are -1)
		mergeQuality = 'Check'

	print '<tr><td><div class="nowrap" style="width:100%%;">%s</div></td><td class="%s">%s</td><td>%s</td><td><a href="%s" target="_blank">Merge %02i</a></td></tr>' % (firstLineOut + '<br>' + alignmentLine + '<br>' + secondLineOut, mergeQuality, mergeQuality, score, actualTEMPFOLDER + 'outFile%02i.txt' % (i+1), i+1)

#	print '<br><a href="%s" target="_blank">Detailed EMBOSS "merger" output for Merge %02i</a> %s (%s)' % (actualTEMPFOLDER + 'outFile%02i.txt' % (i+1), i+1, mergeQuality,score)

print '</table>'

print '</div><br>'

import shutil		# required in block below and later when creating zip archive
if len(referenceSequenceFile) > 0:				# reference sequence specified
	# copy reference sequence file to preserve original
	unalignedFilename = climb.TEMPFOLDER + 'unaligned.fasta'
	shutil.copyfile(refSeqFilename, unalignedFilename)
	theFile = open(unalignedFilename, "a")
	theFile.write('\n')		# ensure at least one newline between reference and merge sequence
	theFile.write('>' + Merge.seq[0]['id'] + '\n')
	theFile.write(Merge.seq[0]['seq'] + '\n')
	theFile.close()

	print '<h3>Reference Alignment</h3>'

	# Specifying a double-length genome as a reference sequence allows for confirmation of merges of full-length sequences
	# which may not merge into the correct orientation because of the specific sequencing primers used for fragments.
	# The final merged sequence will still be in the "wrong" orientation though. This can be corrected by
	# 1. Trimming the aligned sequences and then sliding a specified motif to a specified co-ordinate.
	#    Difficult if the full-length does not include 100% coverage, as then the co-ordinates would not necessarily be accurate.
	# 2. Finding the start of the second copy of the reference sequence in the alignment, and then sliding that residue to position 1.
	#    This could be used even if there was not a 100% match:
	#    1. Store the first n bases in the reference sequence
	#    2. Find the next occurrence of the motif in the alignment
	#    3. Find the corresponding section of the merge sequence which is aligned with this motif
	#    4. Slide that motif to position 1

	print '<a href="#" id="paragraph3Anchor">Click to show/hide reference alignment</a>'
	print '<div id="paragraph3" class="output">'

	alignOut = climb.TEMPFOLDER + 'aligned.fasta'
	# clustal = "clustalw -infile=%s -outfile=%s -align -output=fasta" % (unalignedFilename, alignOut)
	# muscle is faster
	# clustalw does not work if only one reference sequence is provided, as there are only two sequences in total
	# "stable" below preserves input sequence order in output alignment
	#
	# allowing only one reference sequence to be specified (reading the first sequence in the file) means that 'water' can be used (pairwise)
	# and the alignment output of water can be output to the screen, to show similarity between the sequences
	# however, it may be useful to allow more than one reference sequence; for example, deletions may be known to occur in some sequences,
	# so allowing multiple reference sequences means that a few common sequences can be specified.
	# however, this means that 'muscle' or something similar must be used; it may therefore be better to highlight bases with a foreground
	# or background colour, rather than a consensus/similarity symbol. this creates large HTML files with very long lines containing
	# all "span" elements for each base for each sequence; this could be avoided with JavaScript code to do this on the fly, but then the CGI
	# script will have to output the JavaScript code, which creates additional complexity in implementation, testing and debugging.
	# Another alternative would be to allow multiple reference sequences, but to generate the consensus/similarity indicators in the Python script
	# and output these to the output page.
	# This allows multiple sequences to be specified (increased functionality) and avoids large "span" sections or JavaScript.
	muscle = "muscle -in %s -out %s -quiet" % (unalignedFilename, alignOut)

	theFile = os.popen(muscle)
	for l in theFile:
		pass
	theFile.close()

	try:
		theFile = open(alignOut, 'r')
	except IOError:
		htmlPrint('No alignment available.')
	else:
		theFile.close()
		print '<div class="nowrap">'
		# muscle adds newlines to the completed alignment, breaking it into shorter lines
		# easiest way to work around this is to open the file with climb
		out = climb.Sequence()
		out.load(alignOut)

		# determine consensus indicators (mismatches)
		mismatches = ''
		matchCount = 0
		mismatchCount = 0
		for i in range(out.seqLength()[0][1]):		# aligned, so all the same length
			# instead of creating a string of all the bases at each position over all sequences,
			# rather just iterate over all the bases until the first mismatch is found
			# this is a quick and simple mismatch indicator
			j = 0
			mismatchFound = False
			while (j < out.seqCount()-1) and	 (not mismatchFound):
				# print 'i = %i    j = %i    mismatches = %s<br>' % (i, j, mismatches)
				if (out.seq[j]['seq'][i]) != (out.seq[j+1]['seq'][i]):
					mismatchFound = True
				j += 1
			if mismatchFound:
				mismatches += '&nbsp;'		# or '^'
				mismatchCount += 1
			else:
				mismatches += '|'
				matchCount += 1

		widthLimit = 30
		for i in out.seq:
			ttt = i['id'].rstrip()[:widthLimit]	# limit width of FASTA ID for output purposes
			for j in range(widthLimit-len(ttt)):	# bad with _ as leading spaces (%30s) are ignored in html
				ttt += '&nbsp;'
			outLine = '%s %s' % (ttt, i['seq'])
			htmlPrint(outLine)
		ttt = 'Similarity indicator'[:widthLimit]
		for j in range(widthLimit-len(ttt)):
			ttt += '&nbsp;'
		htmlPrint('%s %s' % (ttt, mismatches), 3)
		tempTotal = mismatchCount + matchCount
		htmlPrint('Mismatches: %i/%i (%3.2f%%)' % (mismatchCount, tempTotal, float(mismatchCount)/tempTotal*100))
		htmlPrint('Matches: %i/%i (%3.2f%%)' % (matchCount, tempTotal, float(matchCount)/tempTotal*100))

		# --------------------
		# If double-length reference sequence specified, look in alignment for where the merge
		# sequence aligns with the second instance of the start of the reference sequence
		# and slide the merged sequence accordingly.
		#
		# Decided against this solution as it is heavy and creates difficulties
		# showing the alignment and preparing the files, as the alignment would have to be repeated
		# after the merged seqence has been slid. A better solution is to offer a slide option on the
		# input screen. This means that the merged sequence would be slide and then aligned with a
		# single-length reference sequence. If the merged sequence is over-length, this would have to be
		# trimmed manually; however, such trimming is easier for the end-user than sliding a merged sequence.
		# The length of the merged sequence may not be exactly the same as the expected/reference sequence.
		# Hence, any co-ordinates used may not yield a sequence in the expected position/orientation.
		#
		# for i in out.seq:
		# 	if i['id'] != Merge.seq[0]['id']:	# one of the reference sequences
		#							# reference sequence provided as single- or double-length?
		# slideMotifLength = 25
		# if R.seqLength()[0][1] < slideMotifLength:
		# 	slideMotifLength = R.seqLength()[0][1]		# length of first sequence
		# slideMotif = R.seq[0]['seq'][:slideMotifLength]		# assumes reference sequence is >= 25 bases

		out.unload()

	print '</div>'
	print '</div>'

print '<br>'

# --------------------
# Prepare all files and archive into ZIP for download
import zipfile
zipFolder = 'FMT ' + datetime.datetime.now().strftime("%Y%m%d %Hh%Mm%Ss")
os.chdir(climb.TEMPFOLDER)
os.mkdir(zipFolder)		# avoid making a tar-bomb by placing all files in a timestamped folder
				# files are small, so creating "extra" copies is a better/safer alternative
				# hyperlinks on the output page to the originals will not be broken
os.chdir(zipFolder)
zipList = []
describeList = []
# **1** original sequences: in Merge[1], etc. -- write each one out separately? lastFragment+1

for i in range(1, lastFragment+2):

	if fragABIF[i-1]:		# only output raw / untrimmed for chromatograms
		sourceName = 'Source%02i.fasta' % i
		zipList.append(sourceName)
		describeList.append('%s\tInput sequence in FASTA format' % sourceName)
		f = open(sourceName, 'w')			# use original filename without extension?
		f.write('>' + MergeContainer[i].seq[0]['id'] + '_Source\n')
		f.write(MergeContainer[i].originalCalls + '\n')
		f.close()

	sourceName = 'ToMerge%02i.fasta' % i
	zipList.append(sourceName)
	describeList.append('%s\tSequence used for merge in FASTA format (trimmed if from chromatogram input)' % sourceName)
	f = open(sourceName, 'w')
	f.write('>' + MergeContainer[i].seq[0]['id'] + '_ToMerge\n')
	f.write(MergeContainer[i].seq[0]['seq'] + '\n')
	f.close()
# **2** trimmed sequences: all in tempF[] -- Merge...1...2...
#		filenames stored in method, so re-create with better names above
# **3** merged sequence: Merge0....
fileToMove = os.path.basename(Merge.fileName)
shutil.copy2('../' + fileToMove, '.')		# copy the file from the parent folder
zipList.append(fileToMove)
describeList.append('%s\tFinal merged sequence in FASTA format' % fileToMove)
# **4** merger output files: outFile01...
for i in range(lastFragment):
	fileToMove = 'outFile%02i.txt' % (i+1)
	zipList.append(fileToMove)
	describeList.append('%s\tmerger output' % fileToMove)
	shutil.copy2('../' + fileToMove, '.')		# copy the file from the parent folder

if len(referenceSequenceFile) > 0:
	# **5** reference sequence
	zipList.append(referenceSequenceFile)
	describeList.append('%s\tSupplied reference sequence' % referenceSequenceFile)
	shutil.copy2('../' + referenceSequenceFile, '.')		# copy the file from the parent folder
	# **6*  unaligned.fasta
	zipList.append('unaligned.fasta')
	describeList.append('%s\tFinal merged sequence and reference sequence before alignment' % 'unaligned.fasta')
	shutil.copy2('../unaligned.fasta', '.')		# copy the file from the parent folder
	# **7** aligned.fasta
	zipList.append('aligned.fasta')
	describeList.append('%s\tFinal merged sequence and reference sequence after alignment' % 'aligned.fasta')
	shutil.copy2('../aligned.fasta', '.')		# copy the file from the parent folder

# **8** README
f = open('README.txt','w')
f.write('README\n')
f.write('------\n\n')
f.write('Fragment Merger Tool (HVDR)\n')
f.write('---------------------------\n\n')
f.write('This archive contains all files from the merge\n\n')
f.write('Run completed %s\n\n' % datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
for l in describeList:
	f.write(l + '\n')
f.write('\n-----------\n')
f.close()
zipList.append('README.txt')

# for l in zipList:
# 	print l, '<br>'

for i in range(len(zipList)):
	zipList[i] = zipFolder + '/' + zipList[i]		# cannot just add folder (zipfile), so add each file with path

os.chdir(climb.TEMPFOLDER)
actualZipName = zipFolder + '.zip'
Z = zipfile.ZipFile(actualZipName, 'w')
for zipName in zipList:
	Z.write(zipName)
Z.close()

print '<h3>Archive Download</h3>'
print '<a href="#" id="paragraph4Anchor">Click to show/hide archive download link</a>'
print '<div id="paragraph4" class="output">'

saveFilename = '/tmp/Fragment-Merger-' + randomToken + '/' + actualZipName
htmlPrint('<table cellpadding="5"><tr><td><a href="%s"><div><img src="/res/zip1.png" style="border-style: none; width: 48px; height: 48px"></a></td><td>Right-click icon and select "Save (Link) As" to download a ZIP file archive of all files. The filename is <span class="nowrap">%s</span></td></tr></table>' % (saveFilename, actualZipName))

print '</div><br><br>'

# for i in range(len(fileList)):
# 	filename = "/tmp/Babylon" + randomToken + "/" + fileList[i]
# 	print '<br><a href="%s">Download</a> %s FASTA file<br>' % (filename, proteinList[i])

# filename = "/tmp/Babylon" + randomToken + "/Babylon.zip"
# print '<br><a href="%s">Download</a> all of the above in one ZIP file<br>' % (filename)
# --------------------------------------------------------------------

print '<div class="dotted_red" style="border-width:1.5px; color:darkred; font-size:12px; font-family:monospace; margin-top:5px" align="center">'
print '<b>IT IS ADVISABLE TO CHECK ALL CHROMATOGRAM TRACES CAREFULLY TO ENSURE THAT THE BASE CALLS WHICH HAVE BEEN EXTRACTED ARE ACCURATE AND TO RESOLVE ANY AMBIGUOUS BASE CALLS</b></div>'

# print '<h6><a href="http://commons.wikimedia.org/wiki/File:Merovingian_Bee.svg">Bee</a> &#8226; <a href="http://commons.wikimedia.org/wiki/File:Commons-emblem-success.svg">Green icon</a> &#8226; <a href="http://commons.wikimedia.org/wiki/File:Gnome-emblem-important.svg">Red icon</a> &#8226; <a href="http://commons.wikimedia.org/wiki/File:Commons-emblem-query.svg">Yellow icon</a> &#8226; <a href="http://commons.wikimedia.org/wiki/File:Commons-emblem-notice.svg">Blue icon</a> &#8226; <a href="http://www.alpinemeadow.com/stitchery/weblog/HTML-morsels.html">Percentage bars</a><br></h6>'

print '<h6>Icons from <a href="http://commons.wikimedia.org">Wiki Media Commons</a> &#8226; <a href="http://www.alpinemeadow.com/stitchery/weblog/HTML-morsels.html">Percentage bars</a><br></h6>'

closePage()

# Cleanup /tmp

## MERGE IMAGE: http://commons.wikimedia.org/wiki/File:Symbol_merge_review.svg

## Parse outFile data and report on output page to user
## ZIP all files for download
## Reveal outFile data -- XML Http Request or jQuery

## Detect if duplicate filenames specified (as in error)
## More error checking

## 1 AC   2 CG   3 GT    merges "good" of 1-->3-->2 which is wrong
## Must look at score (negative?), similarity and identity

