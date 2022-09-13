#########################################################################################
#########################################################################################
#
# File: splitVDJ.R
#
# Script to split barcodes from the genetic 'barcode mouse' construct as generated in
# the lab of Ton Schumacher (NKI, NL) in its remaining constant V, D and J elements and
# the modified elements (additions/deletions) in between those constant parts.
#
# Expected input: a data frame with a column named var, which contains the barcode names
# 								(a specific part of the sequence)
#
# Output: 				two data frames are created (add.del.ok and add.del.err), which contain
# 								columns with the remaining constant parts and inserted/deleted parts
#
#
# Author: Joost Beltman (LACDR, Leiden University, NL)
# Date: May 15, 2019
#
#########################################################################################
#########################################################################################


#load required packages
library(seqinr)
library(Biostrings)

#define V, D and J input element
Vpart = "TCCAGTAG"
Vpart.s = s2c(Vpart)
Dfwd = "TCTACTATCGTTACGAC"
Dfwd.s = s2c(Dfwd)
Dinv = "GTCGTAACGATAGTAGA"
Dinv.s = s2c(Dinv)
Jpart = "GTAGCTACTACCG"
Jpart.s = s2c(Jpart)


#make data frame for additions and deletions when compared to input sequence
add.del = as.data.frame( matrix(nrow = nrow(halfsample), ncol = 13) )
names(add.del) = c("V", "insVD", "D", "insDJ", "J", "delV", "delDl", "delDr", "delJ", "D.fwdorinv", "errV", "errD", "errJ")

#calculate barcode totals per mouse
rownames(add.del) = halfsample$var

#split barcodes in individual bases
seqs = rownames(add.del)

if (length(seqs) > 0) {
	for (k in 1:length(seqs)) {
	#for (k in 1:100) {

		#get and split current sequence in list
		curseq = seqs[k]
		curseq.s = s2c(curseq)

		if (length(curseq.s) > 0) { #test if empty sequence

			#find part of V that overlaps	
			overlapV = 0
			for (i in 1:min(length(curseq.s),length(Vpart.s)) ) {

				if (curseq.s[i] == Vpart.s[i]) {
					overlapV = overlapV+1
				} else {
					break
				}
			}

			curseq.V = substr(curseq, 1, overlapV)
		}	else { #empty sequence

			curseq.V = ""

		}

		#remove detected V part from sequence	
		curseq.DJ = substr(curseq, overlapV+1, nchar(curseq))
		curseq.DJ.s = s2c(curseq.DJ)

		if (length(curseq.DJ.s) > 0) {

			#find part of J that overlaps (note that in the loop we now need to start at end)
			overlapJ = 0
			testlen = min(length(Jpart.s), length(curseq.DJ.s))
			for (i in length(Jpart.s):(length(Jpart.s) - testlen + 1)) {

				#find index in sequence to compare to index in J part
				j = length(curseq.DJ.s) - (length(Jpart.s)-i)
				if (curseq.DJ.s[j] == Jpart.s[i]) {
					overlapJ = overlapJ+1
				} else {
					break
				}
			}

			curseq.J = substr(curseq, nchar(curseq)-overlapJ+1, nchar(curseq))

		} else {

			curseq.J = ""

		}

		#remove detected J part from sequence	
		curseq.D = substr(curseq.DJ, 1, nchar(curseq.DJ)-overlapJ)

		if (nchar(curseq.D) > 0) {

			#match first to forward D sequence
			#start with longest possible sequences to find match, and decrease each time by 1
			err.fwd=-1 # to see whether a match was found
			for (l.fwd in min(nchar(curseq.D), nchar(Dfwd)):1) {

				perms.fwd = {}
				for (j in 1:(nchar(Dfwd)-l.fwd+1)) {
					perms.fwd = c(perms.fwd, substr(Dfwd, j, j+l.fwd-1))
				}

				trymatch.fwd = perms.fwd

				#first search exact match	for current length		
				for (i in 1:length(trymatch.fwd)) {

					#find position at which match occurs (-1 for no match)
					matchpos.fwd = regexpr(trymatch.fwd[i], curseq.D, fixed=T)

					#matchpos.fwd = matchPattern(trymatch.fwd[i], curseq.D, max.mismatch=0, with.indels=F)
					#if match found, step out of loop
					if (matchpos.fwd>0) {
						#startpos = start(matchpos.fwd)
						#endpos = end(matchpos.fwd)					
						startpos = matchpos.fwd
						endpos = matchpos.fwd+l.fwd-1

						err.fwd = 0 #remember wether match contains an error
						break
					}
				}

				#if no exact match, allow for one mismatch in middle (last two at both sides must match)
				if (err.fwd<0 & l.fwd >= 5) {

					for (i in 1:length(trymatch.fwd)) {

						matchpos.fwd = matchPattern(trymatch.fwd[i], curseq.D, min.mismatch=1, max.mismatch=1, with.indels=F)

						#if match found, step out of loop
						if (length(matchpos.fwd)>0) {

							startpos = start(matchpos.fwd)
							endpos = end(matchpos.fwd)

							if (substr(curseq.D, startpos, startpos+1) == substr(trymatch.fwd[i], 1, 2) &
									substr(curseq.D, endpos-1, endpos) == substr(trymatch.fwd[i], l.fwd-1, l.fwd)) {

								err.fwd = 1 #remember wether match contains an error
								break

							}
						}

					}
				}

				#if match found, step out of loop
				if (err.fwd>=0) {

					#remember D part in curseq that was matched
					curseq.D.fwdmatch = substr(curseq.D, startpos, endpos)
					curseq.D.fwdmatchVD = substr(curseq.D, 1, startpos-1)
					curseq.D.fwdmatchDJ = substr(curseq.D, endpos+1, nchar(curseq.D))

					#remember parts of D that are left for this match
					Dfwd.leftVD = substr(Dfwd, 1, i-1)
					Dfwd.leftVD.s = s2c(Dfwd.leftVD)
					Dfwd.leftDJ = substr(Dfwd, i+l.fwd, nchar(Dfwd))
					Dfwd.leftDJ.s = s2c(Dfwd.leftDJ)

					break
				}

			}

			#then match to inverted D sequence
			#start with longest possible sequences to find match, and decrease each time by 1
			err.inv=-1 # to see whether a match was found
			for (l.inv in min(nchar(curseq.D), nchar(Dinv)):1) {

				perms.inv = {}
				for (j in 1:(nchar(Dinv)-l.inv+1)) {
					perms.inv = c(perms.inv, substr(Dinv, j, j+l.inv-1))
				}

				trymatch.inv = perms.inv

				#first search exact match	for current length		
				for (i in 1:length(trymatch.inv)) {

					#find position at which match occurs (-1 for no match)
					matchpos.inv = regexpr(trymatch.inv[i], curseq.D, fixed=T)

					#matchpos.inv = matchPattern(trymatch.inv[i], curseq.D, max.mismatch=0, with.indels=F)
					#if match found, step out of loop
					if (matchpos.inv>0) {
						#startpos = start(matchpos.inv)
						#endpos = end(matchpos.inv)					
						startpos = matchpos.inv
						endpos = matchpos.inv+l.inv-1

						err.inv = 0 #remember wether match contains an error
						break
					}
				}

				#if no exact match, allow for one mismatch in middle (last two at both sides must match)
				if (err.inv<0 & l.inv >= 5) {

					for (i in 1:length(trymatch.inv)) {

						matchpos.inv = matchPattern(trymatch.inv[i], curseq.D, min.mismatch=1, max.mismatch=1, with.indels=F)

						#if match found, step out of loop
						if (length(matchpos.inv)>0) {

							startpos = start(matchpos.inv)
							endpos = end(matchpos.inv)

							if (substr(curseq.D, startpos, startpos+1) == substr(trymatch.inv[i], 1, 2) &
									substr(curseq.D, endpos-1, endpos) == substr(trymatch.inv[i], l.inv-1, l.inv)) {

								err.inv = 1 #remember wether match contains an error
								break

							}
						}

					}
				}

				#if match found, step out of loop
				if (err.inv>=0) {

					#remember D part in curseq that was matched
					curseq.D.invmatch = substr(curseq.D, startpos, endpos)
					curseq.D.invmatchVD = substr(curseq.D, 1, startpos-1)
					curseq.D.invmatchDJ = substr(curseq.D, endpos+1, nchar(curseq.D))

					#remember parts of D that are left for this match
					Dinv.leftVD = substr(Dinv, 1, i-1)
					Dinv.leftVD.s = s2c(Dinv.leftVD)
					Dinv.leftDJ = substr(Dinv, i+l.inv, nchar(Dinv))
					Dinv.leftDJ.s = s2c(Dinv.leftDJ)

					break
				}

			}


			#find out whether forward or inverted D fits best when allowing a mismatch at both left and right
			if (l.fwd>4 & ((l.fwd > l.inv) | (l.fwd == l.inv & err.fwd <= err.inv))) { #in case of equal length choose forward match

				curseq.Dnew = curseq.D.fwdmatch

				#find new inserted elements
				curseq.insVD = curseq.D.fwdmatchVD
				curseq.insDJ = curseq.D.fwdmatchDJ

				#find new deleted elements
				curseq.delVD = Dfwd.leftVD
				curseq.delDJ = Dfwd.leftDJ

				#set errors
				add.del[k,]$errD = err.fwd

				#set to forward or unknown
				if ((l.fwd > l.inv) | (l.fwd == l.inv & err.fwd < err.inv)) {
					add.del[k,]$D.fwdorinv = "fwd"
				} else {
					add.del[k,]$D.fwdorinv = "unknown"			
				}
			} else if (l.inv>4) {

				curseq.Dnew = curseq.D.invmatch

				#find new inserted elements
				curseq.insVD = curseq.D.invmatchVD
				curseq.insDJ = curseq.D.invmatchDJ

				#find new deleted elements
				curseq.delVD = Dinv.leftVD
				curseq.delDJ = Dinv.leftDJ

				#set errors
				add.del[k,]$errD = err.inv

				#set to inverted
				add.del[k,]$D.fwdorinv = "inv"	

			} else {

				#too short D fragment matches, so consider D to be empty
				curseq.Dnew = ""

				#find new inserted elements
				if (curseq.V=="" & curseq.J!="") {
					curseq.insVD = substr(curseq.D, 1, nchar(curseq.D))
					curseq.insDJ = ""
				} else if (curseq.J=="" & curseq.V=="") {
					curseq.insVD = ""
					curseq.insDJ = substr(curseq.D, 1, nchar(curseq.D))
				} else {
					curseq.insVD = substr(curseq.D, 1, floor(nchar(curseq.D)/2))
					curseq.insDJ = substr(curseq.D, floor(nchar(curseq.D)/2)+1, nchar(curseq.D))
				}

				#find new deleted elements
				curseq.delVD = substr(Dfwd, 1, nchar(Dfwd)/2)
				curseq.delDJ = substr(Dfwd, nchar(Dfwd)/2+1, nchar(Dfwd))

				#set errors
				add.del[k,]$errD = 0

				#set to inverted
				add.del[k,]$D.fwdorinv = "short"	

			}

			#now the D element is known (including errors allowed at both sides), and the deletions at both sites
			add.del[k,]$D = curseq.Dnew
			add.del[k,]$delDl = curseq.delVD
			add.del[k,]$delDr = curseq.delDJ


			#find out whether V part can be extended if 1 error is assumed
			curseq.VinsVD = paste(curseq.V, curseq.insVD, sep="")
			curseq.VinsVD.s = s2c(curseq.VinsVD)

			#find overlapping part assuming that we start from the V side and move towards D
			#note that first mismatch position is at overlapV+1
			#below we look for a second mismatch position
			maxvalfromV = min(length(Vpart.s), length(curseq.VinsVD.s))
			mismatchposV = overlapV+1
			if ((maxvalfromV - overlapV) > 2) {
				for (overlapVerr in (overlapV+2):maxvalfromV) {
					if (curseq.VinsVD.s[overlapVerr] != Vpart.s[overlapVerr]) {
						mismatchposV = overlapVerr
						break
					}
				}
				#if match occurs until the end of the relevant sequence part, remember this
				if (overlapVerr == maxvalfromV & curseq.VinsVD.s[overlapVerr] == Vpart.s[overlapVerr]) {
					mismatchposV = overlapVerr+1
				}
			}

			#if there is no second mismatch position closeby, we assume there was an error and extend V
			if (mismatchposV > (overlapV+3)) {
				curseq.Vnew = substr(curseq.VinsVD, 1, mismatchposV-1)
				curseq.insVDnewfromV = substr(curseq.VinsVD, mismatchposV, nchar(curseq.VinsVD))
				curseq.delV = substr(Vpart, mismatchposV, nchar(Vpart))
				errorsV=1
			} else { #otherwise no extension
				curseq.Vnew = curseq.V
				curseq.insVDnewfromV = curseq.insVD
				curseq.delV = substr(Vpart, overlapV+1, nchar(Vpart))
				errorsV=0
			}

			#now we know the V part, the deletions from V and the insertions between V and D
			add.del[k,]$V = curseq.Vnew
			add.del[k,]$insVD = curseq.insVDnewfromV
			add.del[k,]$delV = curseq.delV
			add.del[k,]$errV = errorsV

			#find out whether J part can be extended if 1 error is assumed
			curseq.JinsDJ = paste(curseq.insDJ, curseq.J, sep="")
			curseq.JinsDJ.s = s2c(curseq.JinsDJ)

			#find overlapping part assuming that we start from the J side and move towards D
			#note that first mismatch position is at overlapJ+1 (when seen from the back of the sequence)
			#below we look for a second mismatch position
			lenJ = length(Jpart.s)
			lenJinsDJ = length(curseq.JinsDJ.s)
			maxvalfromJ = min(lenJ, lenJinsDJ)
			mismatchposJ = overlapJ+1
			if ((maxvalfromJ - overlapJ) > 2) {
				for (overlapJerr in (overlapJ+2):maxvalfromJ) {
					if (curseq.JinsDJ.s[lenJinsDJ - overlapJerr + 1] != Jpart.s[lenJ - overlapJerr + 1]) {
						mismatchposJ = overlapJerr
						break
					}
				}
				#if match occurs until the end of the relevant sequence part, remember this
				if (overlapJerr == maxvalfromJ & curseq.JinsDJ.s[lenJinsDJ - overlapJerr + 1] == Jpart.s[lenJ - overlapJerr + 1]) {
					mismatchposJ = overlapJerr+1
				}
			}

			#if there is no second mismatch position closeby, we assume there was an error and extend V
			if (mismatchposJ > (overlapJ+3)) {
				curseq.Jnew = substr(curseq.JinsDJ, lenJinsDJ - mismatchposJ + 2, lenJinsDJ)
				curseq.insDJnewfromJ = substr(curseq.JinsDJ, 1, lenJinsDJ - mismatchposJ + 1)
				curseq.delJ = substr(Jpart, 1, lenJ - mismatchposJ + 1)
				errorsJ=1
			} else { #otherwise no extension
				curseq.Jnew = curseq.J
				curseq.insDJnewfromJ = curseq.insDJ
				curseq.delJ = substr(Jpart, 1, lenJ - overlapJ)
				errorsJ=0
			}

			#now we know the J part, the deletions between D and J and the insertions between V and D
			add.del[k,]$J = curseq.Jnew
			add.del[k,]$insDJ = curseq.insDJnewfromJ
			add.del[k,]$delJ = curseq.delJ
			add.del[k,]$errJ = errorsJ

		} else { #empty D element


			add.del[k,]$V = curseq.V
			add.del[k,]$J = curseq.J
			add.del[k,]$D = ""
			add.del[k,]$insVD = ""
			add.del[k,]$insDJ = ""
			add.del[k,]$D.fwdorinv = "short"

			#arbitrary choice of forward element
			add.del[k,]$delDl = substr(Dfwd, 1, nchar(Dfwd)/2)
			add.del[k,]$delDr = substr(Dfwd, nchar(Dfwd)/2+1, nchar(Dfwd))

			#deleted parts of V and J are simply the missing part of the original V and J
			add.del[k,]$delV = substr(Vpart, nchar(curseq.V)+1, nchar(Vpart))
			add.del[k,]$delJ = substr(Jpart, 1, nchar(Jpart)-nchar(curseq.J))

			#no errors at all in this case
			add.del[k,]$errV = 0
			add.del[k,]$errD = 0
			add.del[k,]$errJ = 0

		}
	}
}

add.del.err = add.del[add.del$errV==1 | add.del$errD==1 | add.del$errJ==1 ,]
add.del.ok = add.del[add.del$errV==0 & add.del$errD==0 & add.del$errJ==0 ,]


#write distribution of V, D and J per barcode to file (commented out for usage in combination with tagcleanup script)
#write.table(add.del, paste(data.path.write, "VDJ_indels_perbc.txt", sep=""), sep="\t", quote=F, row.names=T)
