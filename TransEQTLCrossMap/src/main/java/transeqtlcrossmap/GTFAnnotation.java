/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package transeqtlcrossmap;


import umcg.genetica.features.Exon;
import umcg.genetica.features.FeatureType;
import umcg.genetica.features.Gene;
import umcg.genetica.features.Transcript;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;

/**
 * @author Harm-Jan
 */
public class GTFAnnotation extends Annotation {


//	public static void main(String[] args) {
//		String ensemblannotation = "D:\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
//		try {
//			GTFAnnotation g = new GTFAnnotation(ensemblannotation);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//	}


	public GTFAnnotation(String annotation) throws IOException {
		this.annotationLocation = annotation;
		readAnnotation();
	}

	private void readAnnotation() throws IOException {

		System.out.println("Reading annotation: " + annotationLocation);

		TextFile tf = new TextFile(annotationLocation, TextFile.R);
		String ln = tf.readLine();

		strToGene = new HashMap<String, Gene>();

		HashMap<String, Transcript> strToTranscript = new HashMap<String, Transcript>();
		HashMap<String, Exon> strToExon = new HashMap<String, Exon>();

// this all assumes the path is sorted on genomic coordinates..
		HashMap<String, Integer> unknownTypes = new HashMap<>();
		int lnctr = 0;
		while (ln != null) {
			if (ln.startsWith("#")) {
				// skip
			} else {
				GTFLine lineObj = new GTFLine(ln);
				if (lineObj.getType().equals("exon") || lineObj.getType().equals("start_codon") || lineObj.getType().equals("stop_codon") || lineObj.getType().equals("transcript")
						|| lineObj.getType().equals("gene") || lineObj.getType().equals("CDS") || lineObj.getType().equals("UTR")) {
					Gene currentGene = null;
					Transcript currentTranscript = null;

					String geneName = lineObj.getGeneId();
					String transcriptName = lineObj.getTranscriptId();
					Gene tmpGene = new Gene(geneName, lineObj.getChr(), lineObj.getStr());
					if (lineObj.getGeneType() != null) {
						tmpGene.setType(lineObj.getGeneType());
					}
					tmpGene.setGeneSymbol(lineObj.getGeneName());
					if (strToGene.containsKey(tmpGene.getName())) {
						// continue annotating transcripts and exons with this gene
//                System.out.println("getting gene: " + geneName);

						currentGene = strToGene.get(tmpGene.getName());

					} else {
//                System.out.println("new Gene obj");
						currentGene = tmpGene;
						strToGene.put(tmpGene.getName(), currentGene);
					}


					if (!FeatureType.parse(lineObj.getType()).equals(FeatureType.OTHER)) {
						Exon e = new Exon(lineObj.getType(), lineObj.getChr(), lineObj.getStr(), currentGene, lineObj.getStart(), lineObj.getStop());


						e.setName(lineObj.getExonId());

						if (strToTranscript.containsKey(transcriptName)) {
							currentTranscript = strToTranscript.get(transcriptName);
							if (strToExon.containsKey(e.getName())) {
								// no work to be done
								e = strToExon.get(e.getName());
								e.addTranscript(currentTranscript);
								//System.out.println(e.toString());
								//System.out.println("Duplicate entry in path: " + lineObj.toString());
							} else {
								currentTranscript.addExon(e);
								e.addTranscript(currentTranscript);
								strToExon.put(e.getName(), e);
							}

							if (lineObj.getExonNumber() != null) {
								currentTranscript.setExonRank(e, lineObj.getExonNumber());
							}

						} else {
							currentTranscript = new Transcript(lineObj.getTranscriptId(), lineObj.getChr(), lineObj.getStr(), currentGene);
							currentGene.addTranscript(currentTranscript);
							strToTranscript.put(lineObj.getTranscriptId(), currentTranscript);
//                System.out.println(e.toString());
							if (strToExon.containsKey(e.getName())) {
								e = strToExon.get(e.getName());
								e.addTranscript(currentTranscript);
								currentTranscript.addExon(e);
							} else {
								currentTranscript.addExon(e);
								e.addTranscript(currentTranscript);
								strToExon.put(e.getName(), e);
							}
							if (lineObj.getExonNumber() != null) {
								currentTranscript.setExonRank(e, lineObj.getExonNumber());
							}
						}
					}


				} else {
					Integer uobjctr = unknownTypes.get(lineObj.getType());
					if (uobjctr == null) {
						uobjctr = 0;
					}
					uobjctr++;
					unknownTypes.put(lineObj.getType(), uobjctr);
				}
			}
			ln = tf.readLine();
			lnctr++;
			if (lnctr % 1000 == 0) {
				System.out.print(lnctr + " lines parsed\r");
			}
		}
		System.out.print(lnctr + " lines parsed\n");
		tf.close();
		System.out.println(annotationLocation + ", Genes: " + strToGene.size() + "\tTranscripts: " + strToTranscript.size() + "\tExons: " + strToExon.size());

		for (String s : unknownTypes.keySet()) {
			System.out.println("Uknown type found: " + s + ", with occurence: " + unknownTypes.get(s));
		}

		// set the relative start and end positions of each gene and transcript
		genes = strToGene.values();
		genes.forEach(Gene::getBounds);

	}


}
