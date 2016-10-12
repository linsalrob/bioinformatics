/*
Create a fasta file containing all of the hits grouped by things that they hit to
*/

import rob_bio.blast.*;
import rob_bio.*;
import java.io.*;
import java.util.HashMap;

public class CreateFasta {
	public static void main(String args[]) {

		if (args.length < 3) {
			System.err.println("Usage: GroupBlastResults <blast m8 output file> <e-value> <fasta input file> <directory to write to>");
			System.exit(-1);
		}

		GroupBlastResultsByHit hits = new GroupBlastResultsByHit();

		try {
			BufferedReader in = new BufferedReader(new FileReader(args[0]));
			hits = new GroupBlastResultsByHit(in, Double.parseDouble(args[1]), false);
			in.close();
		}
		catch (FileNotFoundException e) {
			System.err.println("File Not Found");
			e.printStackTrace();
		}       
		catch (Exception e) {
			e.printStackTrace();
		}   

		ReadFasta fa = new ReadFasta(args[2]);

		// Create a hash with all the keys and a unique id that we can use to shorten names
		HashMap<String, Integer> idmap = new HashMap();
		try {
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("id.map")));
			int i=0;
			for(String id : fa.ids()) {
				idmap.put(id, ++i);
				out.println(id + "\t" + i);
			}
			out.close();
		}
		catch (Exception e) {
			e.printStackTrace();
		}



		int fileCount = 0;
		for (String k : hits.sortedKeys()) {
			// if (hits.numHits(k) < 2)
			//	continue;

			// System.err.println("Hit " + k + " has  " + hits.numHits(k) + " hits:");
			String filename = args[3] + "/fasta." + ++fileCount;
			try {
				PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
				out.println(">" + idmap.get(k) + "\n" + fa.sequence_of(k));
				for (String match : hits.getHits(k))
					if (k.compareTo(match) != 0)
						out.println(">" + idmap.get(match) + "\n" + fa.sequence_of(match));

				out.flush();
				out.close();
				//System.exit(0);
				}
				catch (SequenceNotFoundException e) {
					System.err.println("BUGGER");
					e.printStackTrace();
				}
				catch (FileNotFoundException e) {
					System.err.println("File Not Found");
					e.printStackTrace();
				}
				catch (Exception e) {
					e.printStackTrace();
				}

			}
		}
	}



