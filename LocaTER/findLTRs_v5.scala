import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.SAMFileWriter
import htsjdk.samtools.SAMFileWriterFactory
import htsjdk.samtools.DefaultSAMRecordFactory
import htsjdk.samtools.ValidationStringency
import htsjdk.samtools.SamReader
import htsjdk.samtools.SAMRecord
import java.io.File
import java.io._
import scala.math._
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet

object findMobileElements{

	def main (args: Array[String]){
		if (args.contains("-h") || args.contains("--help") || args.size < 3){
			println("JAVA_OPTS=-Xmx4g scala -cp htsjdk-1.139.jar:. findMobileElements <SAMPLENAME> <SAMPLEBAM> <REPEATS_DB_UCSC>")
			System.exit(0)
		}

		val Chroms = Array(("chr1",158534110),("chr2",136231102),("chr3",121005158),("chr4",120000601),("chr5",120089316),("chr6",117806340),("chr7",110682743),("chr8",113319770),("chr9",105454467),("chr10",103308737),("chr11",106982474),("chr12",87216183),("chr13",83472345),("chr14",82403003),("chr15",85007780),("chr16",81013979),("chr17",73167244),("chr18",65820629),("chr19",63449741),("chr20",71974595),("chr21",69862954),("chr22",60773035),("chr23",52498615),("chr24",62317253),("chr25",42350435),("chr26",51992305),("chr27",45612108),("chr28",45940150),("chr29",51098607),("chrX",139009144))
		

		val input = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
		//var alignments = input.queryAlignmentStart(Chroms(0)._1,1)
		val inputBreak = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
		//var alignmentsBreak = inputBreak.queryAlignmentStart(Chroms(0)._1,1)
		val repeats = new BufferedReader(new FileReader(new File(args(2))))


		val outFactory = new htsjdk.samtools.SAMFileWriterFactory
		//Write header too new bam
		//val output = outFactory.makeSAMOrBAMWriter(input.getFileHeader,true,new File(s"${args(0)}.bam"))
		//val outFWD = new BufferedWriter(new FileWriter(new File(s"${args(0)}.bedgraph")))
		val outEvent = new BufferedWriter(new FileWriter(new File(s"${args(0)}.event.tab")))
		input.close

		//outFWD.write("track type=bedGraph\n")

		/* Build repeat structures and load data into them*/		

		val repeatsChrom = new HashMap[String,Array[Tuple6[Int,Int,String,String,String,String]]]
		var repList: List[Tuple6[Int,Int,String,String,String,String]] = Nil
		var prevChrom = ""

		while (repeats.ready){
			val tmp = repeats.readLine.split("\t")
			if (prevChrom == tmp(5)){
				repList = (tmp(6).toInt,tmp(7).toInt,tmp(10),tmp(11),tmp(12),tmp(9)) :: repList
			} else {
				repeatsChrom += prevChrom -> repList.reverse.toArray
				repList = (tmp(6).toInt,tmp(7).toInt,tmp(10),tmp(11),tmp(12),tmp(9)) :: Nil
				prevChrom = tmp(5)
			}
		}//end while repeats

		repeatsChrom += prevChrom -> repList.reverse.toArray
		repList = Nil
		prevChrom = ""
		repeats.close

		def isLTR (pos: Int, curSearch: Array[Tuple6[Int,Int,String,String,String,String]]) : Tuple3[Boolean,String,String] = {
			val half = scala.math.floor(curSearch.length/2).toInt
			if (curSearch.length == 0){
					(false,"","")
			} else {
				if (curSearch(half)._1 <= pos && curSearch(half)._2 >= pos) {
					(true,curSearch(half)._6,s"${curSearch(half)._3}-${curSearch(half)._5}")
					} else {
						if (curSearch.length <= 1) {
							(false,"","")
						} else {
						if (curSearch(half)._1 > pos) isLTR(pos,curSearch.slice(0,half)) else isLTR(pos,curSearch.slice(half + 1,curSearch.length))
					}
				}
			}
		} //end is LTR


		/* Begin analysis of Chromosomes for LTRS */


		for (chrs <- Chroms){
		print(" " + chrs._1)
		val input = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
		var alignments = input.query(chrs._1,1,200000000,true)
		var typeEvent = new HashMap[String,Tuple2[Int,Int]]

		var updated = false
		//var LTRpRn, LTRpRp, LTRnRp, LTRnRn = 0
		var RpMp, RpMn, RnMp, RnMn = 0
		var fwdP, revP, fwdS, revS, firstEnd, lastStart = 0
		var hetBridge, readDepth = 0
		var splitEnd = new HashSet[Int]
		var windowDepth = 0
		var allIPP = 0
		var candidateWindowStart, candidateWindowEnd = 0
		var windowBoo = false
		var lastfwd = false
		var c5pSplit = new HashMap[Int,Int]
		var c3pSplit = new HashMap[Int,Int]

		def directionLTR (LTR: String, uniqueNeg: Boolean, mateNeg: Boolean) : Unit = {
			if (uniqueNeg){ 
				// Negative Strand Read
					if (LTR == "+"){
						//LTR is on + strand	
						if (mateNeg) RnMn += 1 else RnMp += 1 
						// Is the mate negative or positive
					} else {
						// LTR is on - strand
						if (mateNeg) RnMp += 1 else RnMn += 1
					}
				} else {
					// Positive Strand Read
					if (LTR == "+"){
						//LTR is on + strand	
						if (mateNeg) RpMn += 1 else RpMp += 1 
						// Is the mate negative or positive
					} else {
						// LTR is on - strand
						if (mateNeg) RpMp += 1 else RpMn += 1
					}
				}
		}

		def updateCounts(breaks: HashSet[Int]) : Unit = {
			if (breaks.size == 2){
				val inputBreak = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
				var breakpoints = breaks.toArray.sorted
				/* Extract reads covering break point if possible */
				//println("Breakpoint refinement")
				val break = inputBreak.queryOverlapping(chrs._1,breakpoints(0),breakpoints(1))
				while (break.hasNext){
					val tmpBreak = break.next
					if ((breakpoints.contains(tmpBreak.getAlignmentStart) || breakpoints.contains(tmpBreak.getAlignmentEnd)) && tmpBreak.getProperPairFlag == true && tmpBreak.getMappingQuality == 60 && tmpBreak.getCigarString.contains("S")) {
						if (tmpBreak.getReadNegativeStrandFlag == false) fwdS += 1 else revS += 1
					}
				}
			}
			inputBreak.close
		}

		def delDup(c5p : HashMap[Int,Int], c3p: HashMap[Int,Int]) : String ={
			if (c5p.size >= 1 && c3p.size >= 1){
					val best5p = c5p.toSeq.sortBy(-_._2).toArray.apply(0)._1
					val best3p = c3p.toSeq.sortBy(-_._2).toArray.apply(0)._1
					if (best5p < best3p) {
						"DUP"
					} else {
						if (best5p > best3p) "DEL" else if (best5p == best3p) "EQL" else "UNK"
					}
				} else{
					"UNK"
				}
		}


		def polyAT(cig: String, read: String): Float = {
			if (cig.indexOf("S") < cig.indexOf("M")) {
				val idx = cig.split("S").apply(0).toInt
				read.substring(0,idx).count(s => s == 'A' || s == 'T')/read.substring(0,idx).length.toFloat
			 } else { 
				val idx =  cig.split("M").apply(0).toInt
				read.substring(idx,read.length).count(s => s == 'A' || s == 'T')/read.substring(idx,read.length).length.toFloat
			}
		}


		def findOthers(breaks: HashSet[Int]) : Unit = {
			val inputBreak = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
			var breakpoints = breaks.toArray.sorted
				/* Extract reads covering break point if possible */
				//println("Breakpoint refinement")
			if (breakpoints.size == 2){
				val break = inputBreak.queryOverlapping(chrs._1,breakpoints(0),breakpoints(1))
				while (break.hasNext){
					val tmpBreak = break.next
					if (tmpBreak.getAlignmentStart < breakpoints(0) && tmpBreak.getAlignmentEnd > breakpoints(1)){
						hetBridge += 1
					}
					readDepth += 1
					if (tmpBreak.getMappingQuality >= 20 && tmpBreak.getMappingQuality <= 60 && ( ( tmpBreak.getCigarString.contains("S") && tmpBreak.getStringAttribute("SA") == null) || tmpBreak.getCigarString.contains("H") ) && repeatsChrom.contains(tmpBreak.getReferenceName) && repeatsChrom.contains(tmpBreak.getMateReferenceName)) {
						if (tmpBreak.getProperPairFlag == true){
							if (tmpBreak.getReadNegativeStrandFlag == false) {
								splitEnd += tmpBreak.getAlignmentStart
								fwdS += 1
								if (c5pSplit.contains(tmpBreak.getAlignmentEnd)) c5pSplit(tmpBreak.getAlignmentEnd) = c5pSplit(tmpBreak.getAlignmentEnd) + 1 else c5pSplit += tmpBreak.getAlignmentEnd -> 1
							} else {
								if (c3pSplit.contains(tmpBreak.getAlignmentStart)) c3pSplit(tmpBreak.getAlignmentStart) = c3pSplit(tmpBreak.getAlignmentStart) + 1 else c3pSplit += tmpBreak.getAlignmentStart -> 1
								splitEnd += tmpBreak.getAlignmentEnd
								revS += 1
							}
						} else {
							val curMateLTRcheck = isLTR(tmpBreak.getMateAlignmentStart,repeatsChrom(tmpBreak.getMateReferenceName))
							val curMainLTRcheck = isLTR(tmpBreak.getAlignmentStart,repeatsChrom(tmpBreak.getReferenceName))
							if (curMainLTRcheck._1 == false && curMateLTRcheck._1 == true){							
								if (tmpBreak.getReadNegativeStrandFlag == false) {
									splitEnd += tmpBreak.getAlignmentEnd
									fwdS += 1
									if (c5pSplit.contains(tmpBreak.getAlignmentEnd)) c5pSplit(tmpBreak.getAlignmentEnd) = c5pSplit(tmpBreak.getAlignmentEnd) + 1 else c5pSplit += tmpBreak.getAlignmentEnd -> 1
								} else {
									if (c3pSplit.contains(tmpBreak.getAlignmentStart)) c3pSplit(tmpBreak.getAlignmentStart) = c3pSplit(tmpBreak.getAlignmentStart) + 1 else c3pSplit += tmpBreak.getAlignmentStart -> 1
									splitEnd += tmpBreak.getAlignmentStart
									revS += 1
								}
							}
						}
					}
				}
			}
			inputBreak.close
		}


		/* Process all alignments on a Chromosome */
		while (alignments.hasNext){
			val tmp = alignments.next		

			/* If the next read is outside of the current region then write out region and advance */
			if (candidateWindowEnd != 0 && tmp.getAlignmentStart >= candidateWindowEnd){
				val targets = splitEnd.toArray.sorted
				if (splitEnd.size == 2 && fwdP > 0 && revP > 0) findOthers(splitEnd) else if (splitEnd.size == 1 && fwdP > 2 && revP > 2) findOthers(HashSet(targets(0) - 10, targets(0) + 10)) else if (splitEnd.size == 0 && fwdP > 4 && revP > 4) findOthers(HashSet(firstEnd,lastStart))
				val splitDif = if (splitEnd.size == 2) scala.math.abs(splitEnd.toArray.apply(0) - splitEnd.toArray.apply(1)) else 0
				if (fwdP >= 3 && revP >= 3 && (fwdS + revS) > 0 && splitDif <= 20) {
					//println(s"Criteria for print ${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}")
					val best5p = if (c5pSplit.size >= 1) c5pSplit.toSeq.sortBy(-_._2).toArray.apply(0)._1 else -1
					val best3p = if (c3pSplit.size >= 1) c3pSplit.toSeq.sortBy(-_._2).toArray.apply(0)._1 else -1
					//if (splitEnd.size == 2 ) updateCounts(splitEnd) else if (splitEnd.size == 1) updateCounts(HashSet(targets(0) - 10, targets(0) + 10))
					//println(s"PASSED ${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}")
					val ornt = s"${RpMp}:${RpMn}:${RnMp}:${RnMn}"
					outEvent.write(s"${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}\t${splitDif}\t${ornt}\t")
					typeEvent.foreach(s => outEvent.write(s"${s._1},${s._2._1},${s._2._2}:"))
					outEvent.write(s"\t${readDepth}:${hetBridge}\t${windowDepth}\t${allIPP - (fwdP+revP)}\t${(fwdP+revP)/(allIPP.toFloat - (fwdP+revP).toFloat + 0.000001)}\t${delDup(c5pSplit,c3pSplit)}\t${best5p}:${best3p}\n")
					//outFWD.write(s"${chrs._1}\t${candidateWindowStart}\t${candidateWindowEnd}\t${fwdP + fwdS + revP + revS}\n")
				} // end of write

				candidateWindowEnd = 0
				candidateWindowStart = 0
				windowBoo = false
				fwdP = 0
				fwdS = 0
				revS = 0
				revP = 0
				allIPP = 0
				splitEnd = new HashSet[Int]
				typeEvent = new HashMap[String,Tuple2[Int,Int]]
				RpMn = 0
				RpMp = 0
				RnMp = 0
				RnMn = 0
				c5pSplit = new HashMap[Int,Int]
				c3pSplit = new HashMap[Int,Int]
				firstEnd = 0
				lastStart = 0
				lastfwd = false				
				readDepth = 0
				hetBridge = 0
				windowDepth = 0
			} // end of window work


			/* PRoper window identification, Check if on good Chromosomes */
			if (repeatsChrom.contains(tmp.getReferenceName) && repeatsChrom.contains(tmp.getMateReferenceName)){
				if (windowBoo) windowDepth += 1
				if (tmp.getProperPairFlag == false && tmp.getMappingQuality >= 20 && tmp.getMappingQuality <= 60 ){
					allIPP += 1
					val curMateLTRcheck = isLTR(tmp.getMateAlignmentStart,repeatsChrom(tmp.getMateReferenceName))
					val curMainLTRcheck = isLTR(tmp.getAlignmentStart,repeatsChrom(tmp.getReferenceName))

					if (curMainLTRcheck._1 == false && curMateLTRcheck._1 == true){

						if (windowBoo == false){
							if (tmp.getReadNegativeStrandFlag == false){
								firstEnd = tmp.getAlignmentEnd
								//output.addAlignment(tmp)
								//println(" create new window " + tmp.getAlignmentStart)
								windowBoo = true
								windowDepth += 1
								fwdP += 1
								lastfwd = true
								directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
								typeEvent += curMateLTRcheck._3 -> (1,0)
								candidateWindowStart = ((tmp.getAlignmentStart / 1000)*1000)
								// 10300 = /1k*1k = 10,000 - 10300 = -300 vs 200 convert to absolute
								if (abs(candidateWindowStart - tmp.getAlignmentStart) >= abs((candidateWindowStart + 500) - tmp.getAlignmentStart)) candidateWindowStart += 500
								candidateWindowEnd = candidateWindowStart + 1500
								//println(candidateWindowStart + " " + candidateWindowEnd)
							} //end if fwd read

						} else {
							/* Window is True, improperly paired and good mapping quality */
							//output.addAlignment(tmp)
							//?? Should this be here? // lastStart = tmp.getAlignmentStart
							//if (! tmp.getCigarString.contains("S")) { //|| tmp.getSupplementaryAlignmentFlag == false){ Not supllemantary alignment
							if (tmp.getStringAttribute("SA") == null ) { 
								if (tmp.getReadNegativeStrandFlag == false){ // Not sup and on +ve strand
									if (lastfwd) firstEnd = tmp.getAlignmentEnd
									if (typeEvent.contains(curMateLTRcheck._3)) typeEvent(curMateLTRcheck._3) = (typeEvent(curMateLTRcheck._3)._1 + 1, typeEvent(curMateLTRcheck._3)._2) else typeEvent += curMateLTRcheck._3 -> (1,0)
									directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
									fwdP += 1
								} else {
									if (lastfwd) {
										lastStart = tmp.getAlignmentStart
										lastfwd = false
									}
									if (typeEvent.contains(curMateLTRcheck._3)) typeEvent(curMateLTRcheck._3) = (typeEvent(curMateLTRcheck._3)._1, typeEvent(curMateLTRcheck._3)._2 + 1) else typeEvent += curMateLTRcheck._3 -> (0,1)
										directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
										revP += 1
								}//end fwd read npp
							} else { // IS Sup Alignment/splitRead
/*
* Add logic that once we start identifing 3' reads we should not count 5' reads, 
* if we see 5' reads we should close the window and reopen a new one.
*/
								//println(" is getSupplementaryAlignmentFlag")
								if (tmp.getReadNegativeStrandFlag == false){
										if (typeEvent.contains(curMateLTRcheck._3)) typeEvent(curMateLTRcheck._3) = (typeEvent(curMateLTRcheck._3)._1 + 1, typeEvent(curMateLTRcheck._3)._2) else typeEvent += curMateLTRcheck._3 -> (1,0)
										fwdS += 1
										fwdP += 1
										if (c5pSplit.contains(tmp.getAlignmentEnd)) c5pSplit(tmp.getAlignmentEnd) = c5pSplit(tmp.getAlignmentEnd) + 1 else c5pSplit += tmp.getAlignmentEnd -> 1
										directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
										splitEnd += tmp.getAlignmentEnd
									} else {
										if (typeEvent.contains(curMateLTRcheck._3))typeEvent(curMateLTRcheck._3) =  (typeEvent(curMateLTRcheck._3)._1, typeEvent(curMateLTRcheck._3)._2 + 1) else typeEvent += curMateLTRcheck._3 -> (0,1)
										revS += 1
										revP += 1
										if (c3pSplit.contains(tmp.getAlignmentStart)) c3pSplit(tmp.getAlignmentStart) = c3pSplit(tmp.getAlignmentStart) + 1 else c3pSplit += tmp.getAlignmentStart -> 1
										directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
										splitEnd += tmp.getAlignmentStart
									} // end sup is neg read
							} // end Sup Else
						} // end if window true
					} // is Repeat structure of mates fine?

				} // end good read mapping and improperly paired

				//if (windowBoo == true) output.addAlignment(tmp)
			}// End of Good Chromosomes check
		} //End of While Alignments
		/* Have completed a chromosome scan now need to clean up*/
		if (windowBoo){
			val targets = splitEnd.toArray.sorted
			if (splitEnd.size == 2 && fwdP > 0 && revP > 0) findOthers(splitEnd) else if (splitEnd.size == 1 && fwdP > 0 && revP > 0) findOthers(HashSet(targets(0) - 10, targets(0) + 10)) else if (splitEnd.size == 0 && fwdP > 4 && revP > 4) findOthers(HashSet(firstEnd,lastStart))
			val splitDif = if (splitEnd.size == 2) scala.math.abs(splitEnd.toArray.apply(0) - splitEnd.toArray.apply(1))+1 else 0
			if (fwdP >= 3 && revP >= 3 && (fwdS + revS) > 0 && splitDif <= 30) {

					val best5p = if (c5pSplit.size >= 1) c5pSplit.toSeq.sortBy(-_._2).toArray.apply(0)._1 else -1
					val best3p = if (c3pSplit.size >= 1) c3pSplit.toSeq.sortBy(-_._2).toArray.apply(0)._1 else -1
					//val targets = splitEnd.toArray.sorted
					//if (splitEnd.size == 2 ) updateCounts(splitEnd) else if (splitEnd.size == 1) updateCounts(HashSet(targets(0) - 10, targets(0) + 10))
					val ornt = s"${RpMp}:${RpMn}:${RnMp}:${RnMn}"
					outEvent.write(s"${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}\t${splitDif}\t${ornt}\t")
					typeEvent.foreach(s => outEvent.write(s"${s._1},${s._2._1},${s._2._2}:"))
					outEvent.write(s"\t${readDepth}:${hetBridge}\t${windowDepth}\t${allIPP - (fwdP+revP)}\t${(fwdP+revP)/(allIPP.toFloat - (fwdP+revP).toFloat + 0.000001)}\t${delDup(c5pSplit,c3pSplit)}\t${best5p}:${best3p}\n")
					//outFWD.write(s"${chrs._1}\t${candidateWindowStart}\t${candidateWindowEnd}\t${fwdP + fwdS + revP + revS}\n")
			}
		}
		input.close

		} // end for chr from Chroms
		outEvent.close
		//outFWD.close
		//output.close
	} // def main

} //object
