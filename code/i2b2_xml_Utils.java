package NEW_NLP;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;

import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import java.io.File;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
public class xml {


	public static String getnote(File XmlFile,String nodename){
		
	
	StringBuilder note= new StringBuilder("");
	try {
		DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
		DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
		Document doc = dBuilder.parse(XmlFile);
		doc.getDocumentElement().normalize();
		NodeList nList = doc.getElementsByTagName(nodename);
		for (int temp = 0; temp < nList.getLength(); temp++) {
			Node nNode = nList.item(temp);
			note.append(nNode.getTextContent());
			
		}
	    } catch (Exception e) {
		e.printStackTrace();
	    }
	return note.toString();
	}
	public static ArrayList<File> walk( String path ) {
		ArrayList<File> filelist = new ArrayList<File>();
        File root = new File( path );
        File[] list = root.listFiles();
        if (list == null) return filelist;
        
        for ( File f : list ) {
            if ( f.isDirectory() ) {
            }
            else {
            	filelist.add(f.getAbsoluteFile());
            }
        }
		return filelist;
    }
  public static String extract_date(String text){
	  StringBuilder date=new StringBuilder();;
	  Pattern stagepattern = Pattern.compile("\\d+-\\d+-\\d+",Pattern.CASE_INSENSITIVE);
	  Matcher matcher = stagepattern.matcher(text);
	  while (matcher.find()) {
	        date.append(matcher.group());
	    }
	    return date.toString();
	}
	  
  
  public static void main(String argv[]) {
	  String folder = "F:/dropbox2gu/Dropbox (Partners HealthCare)/Phenotyping_Protocol/PheCap_Protocol_R1/2014_i2b2_NLP_set1/data";
      ArrayList<File> filelist = walk(folder);
      String note;
      String date;
      for(File f:filelist){
    	  note = getnote(f,"TEXT").trim();
    	  System.out.println(note);
    	  date = extract_date(note.split("\n")[0]);
    	  System.out.println(date);
      }
      

    
  }

}