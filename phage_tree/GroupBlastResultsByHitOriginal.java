
/*
GroupBlastResultsByHit

Rob Edwards October 1st 2008

Read in blast m8 format and group all the results by things that they hit

*/

import java.util.*;
import java.io.*;

public class GroupBlastResults {

	public static void main(String args[]) {

		//HashMap<String, LinkedList> hash = new HashMap();
		HashMap<String, HashSet> hash = new HashMap();
 		HashMap<String, Integer> count 	 = new HashMap();
		HashSet<String> seen 		 = new HashSet();
	
		double cutoff = 1e-50;

		if (args.length < 1) {
			System.err.println("Usage: GroupBlastResults <blast m8 output file>");
			System.exit(-1);
		}

		try {
			BufferedReader in = new BufferedReader(new FileReader(args[0]));

			String line = null;
			while ((line = in.readLine()) != null) {
				String arr[] = line.split("\t", -1);

				double eval = Double.parseDouble(arr[10]);

				if (eval > cutoff)
					continue;

				
				// Uncomment if you want to count only one occurrence
				// if (seen.contains(arr[0]))
				//	continue;

				seen.add(arr[0]);

				if (hash.containsKey(arr[1])) 
					(hash.get(arr[1])).add(arr[0]);
				else {
					//LinkedList<String> ll = new LinkedList();
					HashSet<String> ll = new HashSet();
					ll.add(arr[0]);
					hash.put(arr[1], ll);
				}

				
			}
		}
		catch (FileNotFoundException e) {
			System.err.println("File Not Found");
			e.printStackTrace();
		}
		catch (Exception e) {
			e.printStackTrace();
		}

		/* 
		
			This routine sorts the list based on the value of the keys 

		*/
		

		
		List keys = new ArrayList(hash.keySet());

		for (String k : (String[]) keys.toArray(new String[0]))
			count.put(k, (hash.get(k)).size());

		final HashMap<String, Integer> countComp = count; // This is for the compare() method
		Collections.sort(keys,
				new Comparator(){
				public int compare(Object left, Object right) {
				String lKey = (String) left;
				String rKey = (String) right;

				return ((Integer)countComp.get(rKey).compareTo((Integer)countComp.get(lKey)));
				}
				});

		for (String k : (String[]) keys.toArray(new String[0])) {
			if (count.get(k) <= 1)
				continue;
			//LinkedList<String> ll = hash.get(k);
			HashSet<String> ll = hash.get(k);
			String[] allhits = (String[]) ll.toArray(new String[0]);
			System.out.print(k + "\t" + count.get(k) + "\t");
			for (String h : allhits)
				System.out.print(h + ", ");
			System.out.println();
		}





	}
}
