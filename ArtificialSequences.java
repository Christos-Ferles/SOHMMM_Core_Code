package sohmmm;

import java.util.Random;

class ArtificialSequences extends Alphabet
{

   /*
     Given a probability (0<=p<=1.0) for the occurrence of a symbol retrieve the index of the
     universeOfDiscourse array that corresponds to the symbol that is supposed to be observed
     based on the provided cumulative probability matrix.
   */
   int getLettersIndex(double cumPropMat[], double probability)
   {
      int pointer = 0;
      boolean found = false;

      do
      {
         if(probability>=cumPropMat[pointer]&&probability<cumPropMat[pointer+1])
            found = true;         
         else
            pointer++;
      }
      while(pointer<size&&!found);

      return pointer;
   }

   /*
     Generate a sequence starting with the symbol at position lettex (of the universeOfDiscourse
     array) and continue by utilizing the entries at the respective row of the passed
     cumPropMat[][]. The example found in the main() method is indicative of the use
     of this method.
    */
   char[] getTransitionsSequence(double cumPropMat[][], int lettex, int minLength, int maxLength)
   {
      Random rand = new Random();
      int actualLength = minLength+rand.nextInt(maxLength-minLength+1);
      char sequence[] = new char[actualLength];
      int letterIndex;

      letterIndex = lettex;
      for(int i = 0;i<sequence.length;i++)
      {
         double p = rand.nextDouble();
         letterIndex = getLettersIndex(cumPropMat[letterIndex], p);
         sequence[i] = letters[letterIndex];
      }

      return sequence;
   }

   public static void main(String args[])
   {
      int sizeClusterA = 30;
      int sizeClusterB = 30;
      int sizeClusterC = 30;
      int cardinality = 10;
      int minL = 500;
      int maxL = 500;

      String obSeqs[] = new String[sizeClusterA+sizeClusterB+sizeClusterC];
      String obIds[] = new String[sizeClusterA+sizeClusterB+sizeClusterC];
      ArtificialSequences as = new ArtificialSequences();
      as.setSize(cardinality);

      /* Biased for Transitions */
      double cumPropMat[][] = {
		 {0.000, 0.000, 0.000, 0.500, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
		 {0.000, 0.000, 0.000, 0.500, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
		 {0.000, 0.000, 0.000, 0.000, 0.000, 0.500, 1.000, 0.000, 0.000, 0.000, 0.000},
		 {0.000, 0.000, 0.000, 0.000, 0.000, 0.500, 1.000, 0.000, 0.000, 0.000, 0.000},
		 {0.000, 0.500, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
		 {0.000, 0.500, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
		 {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
		 {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
		 {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
		 {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000}};
      for(int iter = 0;iter<sizeClusterA;iter++)
      {
         obSeqs[iter] =
            new String(as.getTransitionsSequence(cumPropMat, new Random().nextInt(6), minL, maxL));
//                                                           ^^^^^^^^^^^^^^^^^^^^^^^
         obIds[iter] = ">XC->VB->NM-";
      }

      /* Biased for Emissions */
      double prop[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 25.0, 25.0, 25.0, 25.0};
      as.setProbabilityMatrix(prop);
      for(int iter = 0;iter<sizeClusterB;iter++)
      {
         obSeqs[sizeClusterA+iter] = new String(as.getSequence(minL, maxL));
         obIds[sizeClusterA+iter] = "JKLZ";
      }

      as.setUniformProbabilityMatrix();
      for(int iter = 0;iter<sizeClusterC;iter++)
      {
         obSeqs[sizeClusterA+sizeClusterB+iter] = new String(as.getSequence(minL, maxL));
         obIds[sizeClusterA+sizeClusterB+iter] = "MNBVCXZLKJ";
      }

      IOSequences.writeSequences("FirstExperiment(seqs).txt", obSeqs);
      IOSequences.writeSequences("FirstExperiment(ids).txt", obIds);
   }

}
