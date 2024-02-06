package sohmmm;

import java.util.Random;
import java.io.File;

class HiddenMarkovModel
{

   char v[];
   int N;
   int M;
   double a[][];
   double b[][];
   double pi[];
   double w[][];
   double r[][];
   double u[];

   /*
     Creates an array that obeys standard stochastic constraints according
     to the following relaxations:
     All elements sum up to one (NEARLY).
     All elements are non negative (GREATER THAN 1.0E-10).
   */
   double[] standardStochasticConstraintsInitialization(int cardinality)
   {
      double probabilityDistribution[] = new double[cardinality];
      Random rand = new Random();
      double normalizationConstant;
      boolean nonNegative;

      do
      {
         normalizationConstant = 0.0;
         for(int iter = 0;iter<probabilityDistribution.length;iter++)
         {
            probabilityDistribution[iter] = rand.nextDouble();
            normalizationConstant += probabilityDistribution[iter];
         }
         nonNegative = true;
         for(int iter = 0;iter<probabilityDistribution.length;iter++)
         {
            probabilityDistribution[iter] /= normalizationConstant;
            if(probabilityDistribution[iter]<1.0E-10)
               nonNegative = false;
         }
      }
      while(!nonNegative);

      return probabilityDistribution;
   }

   /*
     Displays the sum of all row elements which must be equal to one (NEARLY),
     and the smallest element (numerical) which must be non negative (GREATER
     THAN 1.0E-10).
   */
   void verifyStandardStochasticConstraints(double matrix[])
   {
      double sum = 0.0;
      double smallestElement = Double.POSITIVE_INFINITY;

      for(int iter = 0;iter<matrix.length;iter++)
      {
         sum += matrix[iter];
         if(matrix[iter]<smallestElement)
            smallestElement = matrix[iter];
      }
      System.out.println("TOTAL SUM:: "+sum+"\t\tMINIMUM VALUE:: "+smallestElement);
   }

   void verifyStandardStochasticConstraints(double matrix[][])
   {
      for(int iter = 0;iter<matrix.length;iter++)
         verifyStandardStochasticConstraints(matrix[iter]);
   }

   /*
     Recalculation of A, B and π according to the normalized exponentials
     reparameterization.
   */
   void normalizedExponentialsReparameterization()
   {
      double expSum;

      for(int i = 0;i<a.length;i++)
      {
         expSum = 0.0;
         for(int l = 0;l<w[i].length;l++)
            expSum += Math.exp(w[i][l]);
         for(int j = 0;j<a[i].length;j++)
            a[i][j] = Math.exp(w[i][j])/expSum;
      }

      for(int j = 0;j<b.length;j++)
      {
         expSum = 0.0;
         for(int l = 0;l<r[j].length;l++)
            expSum += Math.exp(r[j][l]);
         for(int t = 0;t<b[j].length;t++)
            b[j][t] = Math.exp(r[j][t])/expSum;
      }

      expSum = 0.0;
      for(int l = 0;l<u.length;l++)
         expSum += Math.exp(u[l]);
      for(int j = 0;j<pi.length;j++)
         pi[j] = Math.exp(u[j])/expSum;
   }

   /*
     Constructor for a Hidden Markov Model. All the observation symbols (case
     sensitive, no blanks) are being retrieved from the observationSymbols
     String, while the cardinality of the state space N is obtained through
     the stateSpaceCardinaltiy variable. A, B, π are initialized according to
     the reparameterization scheme. The elements of the W, R, U matrices are
     being randomly initialized to values from a Gaussian distribution with
     mean 0.0 and standard deviation 1.0.
   */
   HiddenMarkovModel(String observationSymbols, int stateSpaceCardinality)
   {
      v = observationSymbols.toCharArray();
      N = stateSpaceCardinality;
      M = observationSymbols.length();
      a = new double[N][N];
      b = new double[N][M];
      pi = new double[N];
/*
      for(int i = 0;i<a.length;i++)
         a[i] = standardStochasticConstraintsInitialization(a[i].length);
      for(int j = 0;j<b.length;j++)
         b[j] = standardStochasticConstraintsInitialization(b[j].length);
      pi = standardStochasticConstraintsInitialization(pi.length);
*/
      Random randn = new Random();
      w = new double[N][N];
      for(int i = 0;i<w.length;i++)
         for(int j = 0;j<w[i].length;j++)
            w[i][j] = randn.nextGaussian();
      r = new double[N][M];
      for(int j = 0;j<r.length;j++)
         for(int t = 0;t<r[j].length;t++)
            r[j][t] = randn.nextGaussian();
      u = new double[N];
      for(int j = 0;j<u.length;j++)
         u[j] = randn.nextGaussian();
      normalizedExponentialsReparameterization();
   }

   /*
     Constructor for a Hidden Markov Model. All the parameters of the Hidden
     Markov Model are being retrieved from the seven files which are contained
     in the specified folder. It is mandatory that these seven files have been
     created/constructed with the use of function saveParameters().
   */
   HiddenMarkovModel(String folderName)
   {
      String vN[] = IOSequences.readSequences(folderName+File.separator+"vN");
      v = vN[0].toCharArray();
      N = Integer.parseInt(vN[1]);
      M = v.length;
      a = IOFiles.fileToArray(folderName+File.separator+"a", N, N);
      b = IOFiles.fileToArray(folderName+File.separator+"b", N, M);
      pi = IOFiles.fileToArray(folderName+File.separator+"pi", N);
      w = IOFiles.fileToArray(folderName+File.separator+"w", N, N);
      r = IOFiles.fileToArray(folderName+File.separator+"r", N, M);
      u = IOFiles.fileToArray(folderName+File.separator+"u", N);
   }

   /*
     Calculation of the forward variable.
   */
   void forwardVariable(int o[], double alpha[][])
   {
      for(int i = 0;i<alpha[0].length;i++)
         alpha[0][i] = pi[i]*b[i][o[0]];
      for(int t = 1;t<alpha.length;t++)
         for(int j = 0;j<alpha[t].length;j++)
         {
            double sum = 0.0;
            for(int i = 0;i<a.length;i++)
               sum += alpha[t-1][i]*a[i][j];
            alpha[t][j] = sum*b[j][o[t]];
         }
   }

   /*
     Calculation of the scaled forward variable and the scaling coefficient.
   */
   void forwardVariableScalingCoefficient(int o[], double scaledAlpha[][], double sc[])
   {
      for(int t = 0;t<sc.length;t++)
         sc[t] = 0.0;

      for(int i = 0;i<scaledAlpha[0].length;i++)
      {
         scaledAlpha[0][i] = pi[i]*b[i][o[0]];
         sc[0] += scaledAlpha[0][i];
      }
      sc[0] = 1.0/sc[0];
      for(int i = 0;i<scaledAlpha[0].length;i++)
         scaledAlpha[0][i] *= sc[0];

      for(int t = 1;t<scaledAlpha.length;t++)
      {
         for(int i = 0;i<scaledAlpha[t].length;i++)
         {
            double sum = 0.0;
            for(int j = 0;j<a.length;j++)
               sum += scaledAlpha[t-1][j]*a[j][i];
            scaledAlpha[t][i] = sum*b[i][o[t]];
            sc[t] += scaledAlpha[t][i];
         }
         sc[t] = 1.0/sc[t];
         for(int i = 0;i<scaledAlpha[t].length;i++)
            scaledAlpha[t][i] *= sc[t];
      }
   }

   /*
     Calculation of the backward variable.
   */
   void backwardVariable(int o[], double beta[][])
   {
      for(int i = 0;i<beta[beta.length-1].length;i++)
         beta[beta.length-1][i] = 1.0;
      for(int t = beta.length-2;t>=0;t--)
         for(int i = 0;i<beta[t].length;i++)
         {
            beta[t][i] = 0.0;
            for(int j = 0;j<b.length;j++)
               beta[t][i] += a[i][j]*b[j][o[t+1]]*beta[t+1][j];
         }
   }

   /*
     Calculation of the scaled backward variable.
   */
   void backwardVariableScalingCoefficient(int o[], double sc[], double scaledBeta[][])
   {
      for(int i = 0;i<scaledBeta[scaledBeta.length-1].length;i++)
         scaledBeta[scaledBeta.length-1][i] = sc[scaledBeta.length-1];
      for(int t = scaledBeta.length-2;t>=0;t--)
         for(int i = 0;i<scaledBeta[t].length;i++)
         {
            double sum = 0.0;
            for(int j = 0;j<b.length;j++)
               sum += a[i][j]*b[j][o[t+1]]*scaledBeta[t+1][j];
            scaledBeta[t][i] = sc[t]*sum;
         }         
   }

   /*
     Calculation of the likelihood with respect to the forward/backward
     variables.
   */
   double likelihoodComputation(int o[])
   {
      double alpha[][] = new double[o.length][N];

      forwardVariable(o, alpha);
      double p = 0.0;
      for(int j = 0;j<alpha[alpha.length-1].length;j++)
         p += alpha[alpha.length-1][j];
      return p;
   }
       
   /*
     Calculation of the negative log-likelihood with respect to the scaled
     forward/backward variables in order to be kept within machine precision
     bounds.
   */
   double logLikelihoodComputation(int o[])
   {
      double scaledAlpha[][] = new double[o.length][N];
      double sc[] = new double[o.length];

      forwardVariableScalingCoefficient(o, scaledAlpha, sc);
      double logP = 0.0;
      for(int t = 0;t<sc.length;t++)
         logP -= Math.log10(sc[t]);
      return logP;
   }

   /*
     Identical to the polymorphic logLikelihoodComputation method apart from
     the fact that returns the calculated scaledAlpha and sc parameters for
     future (computationally robust) use.
   */
   double logLikelihoodComputation(int o[], double scaledAlpha[][], double sc[])
   {
      forwardVariableScalingCoefficient(o, scaledAlpha, sc);
      double logP = 0.0;
      for(int t = 0;t<sc.length;t++)
         logP -= Math.log10(sc[t]);
      return logP;
   }

   /*
     Searches the (code) index of the observation character according to the
     observation symbols array.
   */
   int getObservationCode(char symbol)
   {
      int code = 0;

      while(v[code]!=symbol)
         code++;
      return code;
   }

   /*
     Transforms an observation sequence array to an array containing the
     corresponding codes from the observation symbols array.
   */
   int[] encodeObservationSequence(char observationSequence[])
   {
      int o[] = new int[observationSequence.length];

      for(int t = 0;t<o.length;t++)
         o[t] = getObservationCode(observationSequence[t]);
      return o;
   }

   /*
     Stores (in an appropriate format) and afterwards retrieves the observation
     sequences and observation identifiers from the respective files.
   */
   Sequences[] loadSequencesIdentifiers(String sequencesFilename,
                                        String identifiersFilename,
                                        String filenamePrefix)
   {
      String observationSequences[];
      String observationIdentifiers[];
      Sequences sequences[];

      observationSequences = IOSequences.loadSequences(sequencesFilename,
                                                       filenamePrefix+".data");
      observationIdentifiers = IOSequences.loadSequences(identifiersFilename,
                                                         filenamePrefix+".info");
      sequences = new Sequences[observationSequences.length];
      for(int d = 0;d<sequences.length;d++)
      {
         sequences[d] = new Sequences();
         sequences[d].o = encodeObservationSequence(observationSequences[d].toCharArray());
         sequences[d].id = observationIdentifiers[d];
      }
      return sequences;
   }

   /*
     Shuffles the observation sequences (alongside with the observation
     identifiers) in order to achieve a random permutation of the employed
     sequences during the learning process.
   */
   static void shuffleSequences(Sequences sequences[])
   {
      Random rand = new Random();
      int index1, index2;
      Sequences swap;

      for(int iter = 0;iter<10*sequences.length;iter++)
      {
         index1 = rand.nextInt(sequences.length);
         index2 = rand.nextInt(sequences.length);
         swap = sequences[index1];
         sequences[index1] = sequences[index2];
         sequences[index2] = swap;
      }
   }

   /*
     Computation of the W parameters.
   */
   void updateW(double hce, int o[], double alpha[][], double beta[][])
   {
      for(int i = 0;i<w.length;i++)
         for(int j = 0;j<w[i].length;j++)
         {
            double sum = 0.0;
            for(int l = 0;l<alpha.length-1;l++)
               sum += alpha[l][i]*(b[j][o[l+1]]*beta[l+1][j]-beta[l][i]);
            w[i][j] += hce*a[i][j]*sum/likelihoodComputation(o);
         }
   }

   /*
     Computation of the W parameters with the use of the scaled variables.
   */
   void scaledUpdateW(double hce, int o[], double sc[], double scaledAlpha[][], double scaledBeta[][])
   {
      for(int i = 0;i<w.length;i++)
         for(int j = 0;j<w[i].length;j++)
         {
            double sum = 0.0;
            for(int l = 0;l<scaledAlpha.length-1;l++)
               sum += scaledAlpha[l][i]*(b[j][o[l+1]]*scaledBeta[l+1][j]-scaledBeta[l][i]/sc[l]);
            w[i][j] += hce*a[i][j]*sum;
         }
   }

   /*
     Computation of the R parameters.
   */
   void updateR(double hce, int o[], double alpha[][], double beta[][])
   {
      for(int j = 0;j<r.length;j++)
         for(int t = 0;t<r[j].length;t++)
         {
            double sum = 0.0;
            for(int l = 0;l<alpha.length;l++)
               if(o[l]==t)
                  sum += (1.0-b[j][t])*alpha[l][j]*beta[l][j];
               else
                  sum -= b[j][t]*alpha[l][j]*beta[l][j];                  
            r[j][t] += hce*sum/likelihoodComputation(o);
         }
   }

   /*
     Computation of the R parameters with the use of the scaled variables.
   */
   void scaledUpdateR(double hce, int o[], double sc[], double scaledAlpha[][], double scaledBeta[][])
   {
      for(int j = 0;j<r.length;j++)
         for(int t = 0;t<r[j].length;t++)
         {
            double sum = 0.0;
            for(int l = 0;l<scaledAlpha.length;l++)
            {
               double conditionalProbability = scaledAlpha[l][j]*scaledBeta[l][j]/sc[l];
               if(o[l]==t)
                  sum += (1.0-b[j][t])*conditionalProbability;
               else
                  sum -= b[j][t]*conditionalProbability;
            }
            r[j][t] += hce*sum;
         }
   }

   /*
     Computation of the U parameters.
   */
   void updateU(double hce, int o[], double beta[][])
   {
      for(int j = 0;j<u.length;j++)
         u[j] += hce*pi[j]*(b[j][o[0]]*beta[0][j]/likelihoodComputation(o)-1.0);
   }

   /*
     Computation of the U parameters with the use of the scaled variables.
   */
   void scaledUpdateU(double hce, int o[], double scaledBeta[][])
   {
      for(int j = 0;j<u.length;j++)
         u[j] += hce*pi[j]*(b[j][o[0]]*scaledBeta[0][j]-1.0);
   }

   /*
     Implementation of the online gradient descent learning algorithm with
     the use of the scaled variables. The parameters needed are the neighborhood
     function's value and the encoded observation sequence under consideration.
   */
   void scaledOnlineGradientDescentLearning(double hce, int o[])
   {
      double scaledAlpha[][] = new double[o.length][N];
      double sc[] = new double[o.length];
      double scaledBeta[][] = new double[o.length][N];

      forwardVariableScalingCoefficient(o, scaledAlpha, sc);
      backwardVariableScalingCoefficient(o, sc, scaledBeta);

      scaledUpdateW(hce, o, sc, scaledAlpha, scaledBeta);
      scaledUpdateR(hce, o, sc, scaledAlpha, scaledBeta);
      scaledUpdateU(hce, o, scaledBeta);

      normalizedExponentialsReparameterization();
/*
      double alpha[][] = new double[o.length][N];
      double beta[][] = new double[o.length][N];
      forwardVariable(o, alpha);
      backwardVariable(o, beta);
      updateW(hce, o, alpha, beta);
      updateR(hce, o, alpha, beta);
      updateU(hce, o, beta);
*/
   }

   /*
     Identical to the polymorphic scaledOnlineGradientDescentLearning method
     apart from the fact that uses the already calculated scaledAlpha and sc
     parameters, leading thus to a computationally robust approach.
   */
   void scaledOnlineGradientDescentLearning(double hce, int o[], double scaledAlpha[][], double sc[])
   {
      double scaledBeta[][] = new double[o.length][N];

      backwardVariableScalingCoefficient(o, sc, scaledBeta);

      scaledUpdateW(hce, o, sc, scaledAlpha, scaledBeta);
      scaledUpdateR(hce, o, sc, scaledAlpha, scaledBeta);
      scaledUpdateU(hce, o, scaledBeta);

      normalizedExponentialsReparameterization();
   }

   /*
     Transform an one dimensional probability matrix to a cumulative
     probability matrix.
   */
   double[] cumulativeProbabilityTransformation(double probabilityMatrix[])
   {
      double cumulativeProbabililtyMatrix[] = new double[probabilityMatrix.length+1];

      cumulativeProbabililtyMatrix[0] = 0.0;
      int iter;
      for(iter = 1;iter<cumulativeProbabililtyMatrix.length-1;iter++)
         cumulativeProbabililtyMatrix[iter] += probabilityMatrix[iter-1]+cumulativeProbabililtyMatrix[iter-1];
      cumulativeProbabililtyMatrix[iter] = 1.0;
      return cumulativeProbabililtyMatrix;
   }

   /*
     Transform a two dimensional (per row) probability matrix to (per row)
     cumulative probability matrix.
   */
   double[][] cumulativeProbabilityTransformation(double probabilityMatrix[][])
   {
      double cumulativeProbabilityMatrix[][] = new double[probabilityMatrix.length][probabilityMatrix[0].length+1];

      for(int iter = 0;iter<cumulativeProbabilityMatrix.length;iter++)
         cumulativeProbabilityMatrix[iter] = cumulativeProbabilityTransformation(probabilityMatrix[iter]);
      return cumulativeProbabilityMatrix;
   }

   /*
     Finds the pointer (index) that represents the probability in the
     cumulative probability array.
   */
   int getIndex(double probability,double cumulativeProbabilityMatrix[])
   {
      int pointer = 0;
      boolean found = false;

      do
      {
         if(probability>=cumulativeProbabilityMatrix[pointer]&&probability<cumulativeProbabilityMatrix[pointer+1])
         {
            found = true;   
         }
         else
            pointer++;
      }
      while(!found);

      return pointer;
   }

   /*
     Generator of encoded observation sequence with desired length.
   */
   int[] generateEncodedObservationSequence(int observationSequenceLength)
   {
      double transitionProbabilities[][] = cumulativeProbabilityTransformation(a);
      double emissionProbabilities[][] = cumulativeProbabilityTransformation(b);
      double initialStateProbabilities[] = cumulativeProbabilityTransformation(pi);
      int o[] = new int[observationSequenceLength];

      Random rand = new Random();
      int state = getIndex(rand.nextDouble(),initialStateProbabilities);
      for(int iter = 0;iter<o.length;iter++)
      {
         o[iter] = getIndex(rand.nextDouble(),emissionProbabilities[state]);
         state = getIndex(rand.nextDouble(),transitionProbabilities[state]);
      }

      return o;
   }

   /*
     Constructor for a Hidden Markov Model used solely by the replicate
     method. The corresponding class variables are only declared without
     being initialized, thus saving unnecessary computations.
   */
   private HiddenMarkovModel()
   {
   }

   /*
     Creates (constructs) and returns a replica (copy) of this object.
   */
   HiddenMarkovModel replicate()
   {
      HiddenMarkovModel hmm = new HiddenMarkovModel();
      hmm.v = new char[M];
      hmm.N = N;
      hmm.M = M;
      hmm.a = new double[N][N];
      hmm.b = new double[N][M];
      hmm.pi = new double[N];
      hmm.w = new double[N][N];
      hmm.r = new double[N][M];
      hmm.u = new double[N];

      System.arraycopy(v, 0, hmm.v, 0, v.length);
      for(int i = 0;i<a.length;i++)
      {
         System.arraycopy(a[i], 0, hmm.a[i], 0, a[i].length);
         System.arraycopy(w[i], 0, hmm.w[i], 0, a[i].length);
         System.arraycopy(b[i], 0, hmm.b[i], 0, b[i].length);
         System.arraycopy(r[i], 0, hmm.r[i], 0, b[i].length);
      }
      System.arraycopy(pi, 0, hmm.pi, 0, pi.length);
      System.arraycopy(u, 0, hmm.u, 0, pi.length);

      return hmm;
   }

   void displayArray(double matrix[])
   {
      for(int iter = 0;iter<matrix.length;iter++)
         System.out.print(matrix[iter]+"\t");
      System.out.println();
   }

   void displayArray(double matrix[][])
   {
      for(int iter = 0;iter<matrix.length;iter++)
         displayArray(matrix[iter]);
      System.out.println();
   }

   /*
     Display all information relavant to the Hidden Markov Model.
   */
   void displayInfo()
   {
      System.out.print("V = {"+v[0]);
      for(int iter = 1;iter<v.length;iter++)
         System.out.print(", "+v[iter]);
      System.out.println("}");
      System.out.println();

      System.out.println("N = "+a[0].length);
      System.out.println();

      System.out.println("M = "+b[0].length);
      System.out.println();

      System.out.println("A =");
      displayArray(a);

      System.out.println("B =");
      displayArray(b);

      System.out.println("\u03C0 =");
      displayArray(pi);
      System.out.println();

      System.out.println("W =");
      displayArray(w);

      System.out.println("R =");
      displayArray(r);

      System.out.println("U ="); 
      displayArray(u);
      System.out.println();
   }

   /*
     Saves/stores all the parameters of the Hidden Markov Model in seven
     files contained in the specified folder.
   */
   void saveParameters(String folderName)
   {
      new File(folderName).mkdir();

      String vN[] = new String[2];
      vN[0] = new String(v);
      vN[1] = String.valueOf(N);
      IOSequences.writeSequences(folderName+File.separator+"vN", vN);

      IOFiles.arrayToFile(a, folderName+File.separator+"a", false);
      IOFiles.arrayToFile(b, folderName+File.separator+"b", false);
      IOFiles.arrayToFile(pi, folderName+File.separator+"pi", false);

      IOFiles.arrayToFile(w, folderName+File.separator+"w", false);
      IOFiles.arrayToFile(r, folderName+File.separator+"r", false);
      IOFiles.arrayToFile(u, folderName+File.separator+"u", false);
   }

}
