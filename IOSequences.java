package sohmmm;

import java.io.*;
import java.util.Vector;

class IOSequences
{

   /*
      Reads sequences delimited by '\n', the file mustn't contain empty
      spaces or empty lines (consequently the last sequence should end
      with '\n' <EOF>). It is also assumed that the file contains
      at least one sequence.
   */
   static String[] readSequences(String filename)
   {
      Vector<String> stringArray = new Vector<String>(1,1);

      try
      {
         FileReader fr = new FileReader(filename);
         boolean eof = false;
         String str;
         int ascii, index = 0;

         while(!eof)
         {
            str = "";
            while((ascii = fr.read())!='\n'&&ascii!=-1)
               str += Character.toString((char)ascii);
            stringArray.add(index,str);
            index++;
            if(ascii==-1)
               eof = true;
         }
         fr.close();
      }
      catch(IOException e)
      {
         System.out.println("READ ERROR: "+e.toString());
      }

      String sequences[] = new String[stringArray.size() - 1];
      for(int i = 0;i<sequences.length;i++) {
         sequences[i] = (String)stringArray.get(i);
      }
      return sequences;
   }

   /*
      Writes sequences delimited by '\n', it doesn't write empty spaces or
      empty lines. As expected the last sequence ends with '\n' <EOF>.
      It is assumed that at least one sequence is going to be written.
   */
   static void writeSequences(String filename,String sequences[])
   {
      try
      {
         FileWriter fw = new FileWriter(filename);
         int i;

         for(i = 0;i<sequences.length;i++)
         {
            fw.write(sequences[i],0,sequences[i].length());
            fw.write('\n');
         }
         fw.close();
      }
      catch(IOException e)
      {
         System.out.println("WRITE ERROR: "+e.toString());
      }
   }

   /*
     Removes '\r' characters and empty lines from text files. Consequently,
     the last line ends with '\n' <EOF>. Window's text editors use the
     combination '\r' and '\n' to indicate a new line.
   */
   static void standardFormat(String sourceFile, String destinationFile)
   {
      try
      {
         FileReader fr = new FileReader(sourceFile);
         FileWriter fw = new FileWriter(destinationFile);
         boolean eof = false;
         String str;
         int ascii;

         while(!eof)
         {
            str = "";
            while((ascii = fr.read())!='\n' && ascii!=-1)
            {
//               System.out.println(ascii);
               str += Character.toString((char)ascii);
            }
            if(str.compareTo("\r")!=0 && str.compareTo("")!=0)
            {
            	if(str.contains("\r")) {
            		str = str.replace('\r', '\n');
            	}
            	else if(ascii!=-1) {
            		str += Character.toString('\n');
            	}
            	fw.write(str,0,str.length());
            }
            if(ascii==-1)
               eof = true;
         }
         fw.close();
         fr.close();
      }
      catch(IOException e)
      {
         System.out.println("FORMAT ERROR: "+e.toString());
      }
   }   
   
   static void stats(String checkFile)
   {
      try
      {
         FileReader fr = new FileReader(checkFile);
         int ascii;
    	 int rCount = 0, nCount = 0; 

         while((ascii = fr.read())!=-1) {
        	 if (ascii == '\n')
        		 nCount++;
        	 else if (ascii == '\r')
        			 rCount++;
         }
         fr.close();
         
         System.out.println("\\r = " + rCount);
         System.out.println("\\n = " + nCount);         
      }
      catch(IOException e)
      {
         System.out.println("FORMAT ERROR: "+e.toString());
      }
   } 

   /*
     Combines the functionalities of the standardFormat and readSequences
     static methods.
   */
   static String[] loadSequences(String sequenceFile, String formattedFile)                                                                           
   {
      String sequences[];

      standardFormat(sequenceFile, formattedFile);
      sequences = readSequences(formattedFile);
      return sequences;
   }   

}
