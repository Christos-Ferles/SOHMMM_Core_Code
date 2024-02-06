package sohmmm;

import java.util.Random;

/*
Class for generating sequences of symbols (the current options include all uppercase letter
and digit characters but this is something easily expandable). Generating sequences and
estimating related statistics is pretty straightforward and can be accomplished after setting
the size of the "alphabet" and the corresponding observation probability of each symbol.
*/
class Alphabet
{

   final char universeOfDiscourse[] = {
                                       'M', 'N', 'B', 'V', 'C', 'X',
                                       'Z', 'L', 'K', 'J', 'H', 'G',
                                       'F', 'D', 'S', 'A', 'P', 'O',
                                       'I', 'U', 'Y', 'T', 'R', 'E',
                                       'W', 'Q', '0', '9', '8', '7',
                                       '6', '5', '4', '3', '2', '1'
                                      };
   int size;
   char letters[];
   double probabilityMatrix[];

   /*
      Setting the size of the utilized alphabet/universe of discourse.
    */
   boolean setSize(int size)
   {
      if(size>0&&size<=universeOfDiscourse.length)
      {
         this.size = size;
         letters = new char[this.size];
         for(int i = 0;i<letters.length;i++)
            letters[i] = universeOfDiscourse[i];
         return true;
      }
      else
         return false;
   }
   
   /*
      Setting the probability of occurrence for each respective symbol from the employed
      universe of discourse.
    */
   void setProbabilityMatrix(double percentages[])
   {
      probabilityMatrix = new double[letters.length+1];

      probabilityMatrix[0] = 0.0;
      probabilityMatrix[probabilityMatrix.length-1] = 1.0;
      for(int i = 1;i<probabilityMatrix.length-1;i++)
         probabilityMatrix[i] = probabilityMatrix[i-1]+percentages[i-1]/100.0;
//      for(int i = 0;i<probabilityMatrix.length;i++)
//         System.out.println(probabilityMatrix[i]);
   }
   
   /*
      Setting equal occurrence/observation probabilities for each utilized symbol.
    */
   void setUniformProbabilityMatrix()
   {
      double step = 100.0/(double)letters.length;
      double percentages[] = new double[letters.length];
      for(int i = 0;i<percentages.length;i++)
         percentages[i] = step;
      setProbabilityMatrix(percentages);
   }
   
   /*
     Given a probability (0<=p<=1.0) for the occurrence of a symbol retrieve the index of the
     universeOfDiscourse array that corresponds to the symbol that is supposed to be observed
     with this frequency/probability.
    */
   int getLettersIndex(double probability)
   {
      int pointer = 0;
      boolean found = false;

      do
      {
         if(probability>=probabilityMatrix[pointer]&&probability<probabilityMatrix[pointer+1])
         {
//            System.out.println(probability+" => "+letters[pointer]);
            found = true;
         }
         else
            pointer++;
      }
      while(pointer<size&&!found);

      return pointer;
   }
   
   /*
     Generate a sequence based on the designated statistics with length in the
     defined range.
    */
   char[] getSequence(int minLength, int maxLength)
   {
      Random rand = new Random();
      char sequence[];

      int actualLength = minLength+rand.nextInt(maxLength-minLength+1);
      sequence = new char[actualLength];
      for(int i = 0;i<sequence.length;i++)
      {
         double probability = rand.nextDouble();
         sequence[i] = letters[getLettersIndex(probability)];         
      }
      return sequence;
   }
   
   /*
     Generate a sequence based on the designated statistics with length in the
     defined range, without allowing consecutive (double, triple etc.) identical
     symbols.
   */
   char[] getRepetitionlessSequence(int minLength, int maxLength)
   {
      char sequence[];
      double probability;
      Random rand = new Random();
      char letter;

      int actualLength = minLength+rand.nextInt(maxLength-minLength+1);
      sequence = new char[actualLength];
      probability = rand.nextDouble();
      sequence[0] = letters[getLettersIndex(probability)];               
      for(int i = 1;i<sequence.length;i++)
      {
         probability = rand.nextDouble();
         letter = letters[getLettersIndex(probability)];
         if(letter!=sequence[i-1])
            sequence[i] = letter;
         else
            i--;  
      }
      return sequence;
   }
   
   /*
     Given an already created sequence, generate a new one by carrying out changes
     (replacements, insertions, deletions) according to the provided individual
     probabilities.
   */   
   char[] sequenceVariations(char initialSequence[], double replace, double insert, double delete)
   {
      Random rand = new Random();
      String midSequence = "";

      replace /= 100.0;
      insert /= 100.0;
      delete /= 100.0;
      for(int i = 0;i<initialSequence.length;i++)
      {
         double probability = rand.nextDouble();
         if(probability>=0&&probability<replace)
         {
            double pi8anothta = rand.nextDouble();
            int pointer = getLettersIndex(pi8anothta);
            midSequence += String.valueOf(letters[pointer]);
//            System.out.println("REPLACE ("+i+"): "+initialSequence[i]+"->"+letters[pointer]);
         }
         else if(probability>=replace&&probability<replace+insert)
         {
            double pi8anothta = rand.nextDouble();
            int pointer = getLettersIndex(pi8anothta);
            midSequence += String.valueOf(letters[pointer]);
//            System.out.println("INSERT ("+i+"): "+letters[pointer]);
            i--;
         }
         else if(probability>=replace+insert&&probability<replace+insert+delete)
         {
//            System.out.println("DELETE ("+i+"): "+initialSequence[i]);
         }
         else
         {
            midSequence += String.valueOf(initialSequence[i]);
//            System.out.println("COPY ("+i+"): "+initialSequence[i]);
         }     
      }     
      char sequence[] = new char[midSequence.length()];
      for(int i = 0;i<midSequence.length();i++)
         sequence[i] = midSequence.charAt(i);
      return sequence;
   }      
   
   /*
     Retrieve the index of the universeOfDiscourse array that corresponds
     to the symbol under consideration.
   */
   int getIndex(char letter)
   {
      int index = -1;

      for(int i = 0;i<letters.length;i++)
         if(letters[i]==letter)
         {
            index = i;
            break;
         }
      return index;
   }
   
   /*
     Show the sequence with each symbol seperated by an empty space.
    */
   void displaySequence(char sequence[])
   {
      for(int i = 0; i<sequence.length - 1; i++)
         System.out.print(sequence[i] + " ");
      System.out.println(sequence[sequence.length - 1]);
   }
   
   /*
     Show the occurrence frequency for each contained symbol in the sequence.
    */
   void showStatistics(char sequence[])
   {
      int percentages[] = new int[size];

      for(int i = 0;i<percentages.length;i++)
         percentages[i] = 0;
      for(int i = 0;i<sequence.length;i++)
         percentages[getIndex(sequence[i])]++;
      for(int i = 0;i<percentages.length;i++)
         System.out.println(letters[i]+":: "+percentages[i]*100.0/sequence.length+" (%)");
   }

}
