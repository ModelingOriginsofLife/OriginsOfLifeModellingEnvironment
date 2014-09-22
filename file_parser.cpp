#include "includes.h"

unordered_map<string,int> sectionHeaders({
	{ "LIBRARY", SECTION_LIBRARY },
	{ "CONJUGATES", SECTION_CONJUGATES },
	{ "RULES", SECTION_RULES },
	{ "INITIAL", SECTION_INITIAL },
	{ "CONTROL", SECTION_CONTROL },
	{ "BOUNDARY", SECTION_BOUNDARY },
	{ "GEOMETRY", SECTION_GEOMETRY },
	{ "END", SECTION_END },
	{ "CONFIG", SECTION_CONFIG }
});

void parseError(const char *fmt, ...)
{
	va_list argp;
	
	va_start(argp, fmt);		
	vprintf(fmt, argp);
	exit(0);
}

void getCollapsedLine(ifstream& file, string& str)
{
	int comment=0;
	int wscount = 0;
	
	do
	{
		getline(file,str);
		
		if (file.rdstate() & ifstream::badbit) parseError("Unexpected end of file.\n");
		
		wscount = 1; // start with whitespace
		comment = 0;
		for (int i=0;i<str.length();i++)
		{			
			if (comment)
			{
				str.erase(i,1);
				i--;
			}
			else			
			if (str[i]=='#')
			{
				str.erase(i,1);
				i--;
				comment=1;
			}
			else
			if ((str[i]==' ')||(str[i]=='\t'))
			{
				wscount++;
				
				if (wscount>1)
				{
					str.erase(i,1);
					i--;
				}
				else
				{
					str[i] = ' '; // Standardize whitespace					
				}
			} else wscount = 0;
		}
	} while (str.length()==0); // Skip empty lines or lines that are only comment
	
//	printf("%s\n",str.c_str());
}

string getFirstToken(string input, string& token)
{
	int pos = input.find_first_of(" ");
	
	if (pos != string::npos)
	{
		token=input.substr(0,pos);	
		return input.substr(pos+1,input.length()-pos-1);
	}
	else
	{
		token=input;
		
		string empty="";
		return empty;
	}
}

int getSectionID(string str)
{
	if (sectionHeaders.count(str)) return sectionHeaders[str];
	else return -1;
}

void ChemistryComputation::parseLibrary(ifstream& file)
{
	string str;
	int istop = 0;
	
	do
	{
		getCollapsedLine(file,str);
		if (str=="ENDSECTION") istop=1;
		else {
			if (str.length()==1)
			{
				L.atomTypes.push_back(str[0]);
			}
			else
			{
				// label: a,b,c,...
				char label = str[0];
				int pos = str.find_first_of(" ,");
				vector<char> members;
						
				while (pos != string::npos)
				{
					if (pos<str.length()-1)
					{
						if ((str[pos+1]!=' ')&&(str[pos+1]!=','))
						{
							members.push_back(str[pos+1]);
						}
					}
					
					pos = str.find_first_of(" ,",pos+1);
				}
				
				L.atomClasses[label] = members;
			}
		}
	} while (!istop);
}

void ChemistryComputation::parseSingleConjugateSet(ifstream& file, int setidx)
{
	string str, substr;
	int istop = 0;
	
	do
	{
		getCollapsedLine(file,str);
		if (str=="ENDSET") istop=1;
		else
		{
			// src sink
			char label = str[0];
			int pos = str.find_first_of(" ,");
			char target = str[pos+1];
			
			L.conjugates[setidx][(unsigned char)label] = (unsigned char)target;
			
		}
	} while (!istop);
}

void ChemistryComputation::parseConjugates(ifstream& file)
{
	string str, substr;
	int istop = 0;
	int setidx;
	
	do
	{
		getCollapsedLine(file,str);
		str = getFirstToken(str, substr);
		if (substr=="ENDSECTION") istop=1;
		else if (substr == "SET")
		{
			setidx = atoi(str.c_str());
			
			L.expandConjugates(setidx);
			
			parseSingleConjugateSet(file, setidx);
		}
	} while (!istop);
}

void ChemistryComputation::parseSingleRule(ifstream& file)
{
	ReactionRule R;
	string str, substr;
	int istop;
	
	R.rule = "";
	R.k = 1; // Default reaction constant			
	
	do
	{
		getCollapsedLine(file,str);
		str = getFirstToken(str, substr);
		
		if (substr=="RATE")
		{
			R.k = atof(str.c_str());
		} else if (substr=="RULE")
		{
			R.rule = str;
		} else if (substr=="FUNCTION")
		{
		} else if (substr=="CONDITION")
		{
		} else if (substr=="ENDRULE") istop=1;		
	} while (!istop);
			
	reactions.push_back(R);
}

int ChemistryComputation::lookupRegion(string str)
{
	if (!simTemplate.regionMap.count(str))
	{
		simTemplate.regionMap[str] = simTemplate.regions.size();
		
		simTemplate.regions.emplace(simTemplate.regions.end());
		
		return simTemplate.regions.size()-1;
	} else return simTemplate.regionMap[str];
}

void ChemistryComputation::parseInitial(ifstream& file)
{
	int istop = 0;
	string str, substr;
	int place = 0;
	
	do
	{
		getCollapsedLine(file,str);				
		str = getFirstToken(str, substr);
		if (substr=="ENDSECTION") istop=1;
		else if (substr=="LOCATION") { place = lookupRegion(str); }
		else if (substr=="BATH") 
		{
			str = getFirstToken(str, substr); // Compound, number
			int count = atoi(str.c_str());
			
			simTemplate.regions[place].bath.addNewLeaf(substr, count);
		}
		else if (substr=="COMPOUND") 
		{
			str = getFirstToken(str, substr); // Compound, number
			int count = atoi(str.c_str());
			
			simTemplate.regions[place].addCompound(substr, count);			
		}		
	} while (!istop);
}

void ChemistryComputation::parseControl(ifstream& file)
{
	int istop = 0;
	string str, substr;
	int place = 0;
	
	do
	{
		getCollapsedLine(file,str);				
		str = getFirstToken(str, substr);

		if (substr == "ANALYSIS")
		{
			for (int i=0;i<registeredAnalyses.size();i++)
			{
				if (str == registeredAnalyses[i]->analysisType)
				{
					AnalysisRequest *R = registeredAnalyses[i]->clone();
					R->parse(file);
					analyses.push_back(R);
				}
			}
		}
		else if (substr == "SIMULATE")
		{
			for (int i=0;i<registeredSimulations.size();i++)
			{
				if (str == registeredSimulations[i]->simType)
				{
					SimulationRequest *R = registeredSimulations[i]->clone();					
					R->parse(file);
					jobs.push_back(R);
				}
			}
		}
		else if (substr == "LOAD") // Load initial conditions from a file
		{
		}
		else if (substr == "SAVE") // Save system state to a file during the simulation
		{
		}
		else if (substr == "PERTURB") // Alter system state at particular times during simulation runs
		{
		}
		else if (substr == "ENDSECTION") istop=1;
	} while (!istop);
}

void ChemistryComputation::parseConfig(ifstream& file)
{
	int istop = 0;
	string str, substr;
	int place = 0;
	
	do
	{
		getCollapsedLine(file,str);				
		str = getFirstToken(str, substr);
	
		if (substr == "SINGLETRATE")
		{
			if (str == "LENGTH") singletRate = WEIGHT_HEAVYLENGTH;
			else if (str == "MOLECULE") singletRate = WEIGHT_HEAVY;
		}
		else if (substr=="ENDSECTION") istop=1;
	} while (!istop);
}

void ChemistryComputation::parseRules(ifstream& file)
{
	int istop = 0;
	string str;
	
	do
	{
		getCollapsedLine(file,str);				
		if (str=="ENDSECTION") istop=1;
		else if (str=="ADDRULE")
		{
			parseSingleRule(file);
		}
	} while (!istop);
}

ChemistryComputation parseConfigFile(ifstream& file)
{
	ChemistryComputation C;
	string str;
	int stop=0;
	
	do
	{		
		getCollapsedLine(file, str);
		int section = getSectionID(str);

		switch (section)
		{
			case SECTION_LIBRARY:
				C.parseLibrary(file);
				break;
			
			case SECTION_CONJUGATES:
				C.parseConjugates(file);
				break;
			
			case SECTION_RULES:
				C.parseRules(file);
				break;
				
			case SECTION_INITIAL:
				C.parseInitial(file);
				break;
				
			case SECTION_BOUNDARY:
				break;
				
			case SECTION_END:
				stop = 1;
				break;
				
			case SECTION_CONTROL:
				C.parseControl(file);
				break;
				
			case SECTION_CONFIG:
				C.parseConfig(file);
				break;
				
			default:
				parseError("Invalid section header: %s\n", str.c_str());
				break;
		}
	} while (!stop);
	
	for (int i=0;i<C.reactions.size();i++)
	{
		if (C.reactions[i].rule.length())
		{
			C.reactions[i].parseRule(C.L);
		}
	}
	
	return C;
}
