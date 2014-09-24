#include "includes.h"

void Symbol::clear()
{
	isWild = isSingleton = isReversed = isBound = false;
	doConjugate = 0;
	label=' ';	
}

string reverseString(string S)
{
	string R = "";
	
	for (int i=0;i<S.length();i++)
	{
		R = R + S[S.length()-i-1];
	}
	
	return R;
}

vector<vector<string>> SearchSubtree::generateProducts(vector<vector<Symbol>> rules, Library L)
{
	vector<vector<string>> pList;
	
	if (subTrees.size() == 0)
	{
		if (isEnd)
		{
			vector<string> subList;
			
			for (int i=0;i<rules.size();i++)
				subList.push_back(getProduct(rules[i], L));
				
			pList.push_back(subList);
			return pList;
		}
		else
		{
			return pList;
		}
	}
	else
	{		
		for (int i=0;i<subTrees.size();i++) // Descend into the truee
		{
			vector<vector<string>> subList = subTrees[i].generateProducts(rules, L);
			pList.insert(pList.end(), subList.begin(), subList.end());
		}
	}
	
	return pList;
}

string SearchSubtree::getProduct(vector<Symbol> rule, Library L)
{
	int i;
	string productString = "";
	
	for (i=0;i<rule.size();i++)
	{
		if (rule[i].isBound)
		{
			if (!rule[i].isReversed)
				productString = productString + L.conjugateString(bound[rule[i].label], rule[i].doConjugate);
			else
				productString = productString + reverseString(L.conjugateString(bound[rule[i].label], rule[i].doConjugate));
		}
		else
		{
			if (rule[i].isWild)
			{
				// This is a problem! There's no way to know what this wildcard should render as!
				// So lets just do nothing...
			}
			else
			{
				productString = productString + L.applyConjugation(rule[i].label, rule[i].doConjugate);
			}
		}
	}
	
	return productString;
}

void SearchSubtree::parseOnLeaves(vector<Symbol> rule, string matchstr, Library L)
{
	if (subTrees.size() == 0)
	{
		if (isEnd)
		{
			isEnd = 0;
			rulepos = 0; stringpos = 0;
			applyRule(rule, matchstr, L);
		}
	}
	else
	{
		for (int i=0;i<subTrees.size();i++) // Descend into the truee
			subTrees[i].parseOnLeaves(rule, matchstr, L);
	}
}

void Library::expandConjugates(int idx)
{
	if (conjugates.size() > idx) return;
	
	int i;
	
	for (i=conjugates.size();i<=idx;i++)
	{
		vector<char> newConjugate;
		
		newConjugate.resize(256,0);
		
		conjugates.push_back(newConjugate);
	}
}

char Library::applyConjugation(char C, int idx)
{
	if (idx==0) return C;
	
	return (conjugates[idx])[(unsigned char)C];
}

string Library::conjugateString(string S, int idx)
{
	if (idx==0) return S;
	else
	{
		string R = "";
		
		for (int i=0;i<S.length();i++)
		{
			R = R + applyConjugation(S[i],idx);
		}
		return R;
	}	
}

bool Library::isMemberOfSet(char C, char label, int idx)
{
	if (!atomClasses.count(label))
	{
		return (C==applyConjugation(label,idx));
	}	
	
	if (idx==0)
	{
		if (find(atomClasses[label].begin(), atomClasses[label].end(), C)!=atomClasses[label].end())
			return true;
		return false;
	}
	else
	{
		vector<char>::iterator begin, end;
		
		begin=atomClasses[label].begin();
		end=atomClasses[label].end();
		
		for (;begin!=end;++begin)
		{
			if (applyConjugation(*begin, idx) == C)
			{
				return true;
			}
		}
	}
	
	return false;
}

void SearchSubtree::writeValidLeaves()
{
	if (isNull) return;

	if (isEnd)
	{
		printf("Leaf: ");
		unordered_map<char,string>::iterator it;
		
		for (it=bound.begin();it!=bound.end();it++)
			printf("%c: %s, ",it->first, it->second.c_str());
		printf("\n");
		
		return;
	}
	
	for (int i=0;i<subTrees.size();i++)
		subTrees[i].writeValidLeaves();
}

bool SearchSubtree::isMatch(Symbol S, char C, Library L)
{	
	if (S.isSingleton)
	{
		return (C==L.applyConjugation(S.label, S.doConjugate));
	}
	
	if (S.isBound)
	{
		if (bound.count(S.label))
		{
			return (L.applyConjugation(bound[S.label][0], S.doConjugate) == C);
		}
		else
		{
			string str = ""; str = str + C;
			bound[S.label] = str; 
			return true;
		}
	}
	
	return true;
}

void SearchSubtree::applyRule(vector<Symbol> rule, string matchstr, Library L)
{
	int offset=rulepos, newoffset, strpos = stringpos;
	int stop = 0;
	int i;
	
	if (offset>=rule.size())
	{
		if (strpos < matchstr.length())
		{
			isNull = true; // Did not finish matching compound
		}
		else
			isEnd = true;
		return;
	}
	
	do
	{
		if (strpos >= matchstr.length()) 
		{
			if ((rule[offset].isWild)&&(!bound.count(rule[offset].label))) // This is to check if we can skip this character if we're out of string
			{
				branchTree(offset, strpos, rule, matchstr, L);
				for (i=0;i<subTrees.size();i++)
				{
					subTrees[i].applyRule(rule,matchstr,L);
				}
				return;
			} 
			else 
				stop=1; // Exceeded length of target compound
		}
		else
		{
			if (rule[offset].isWild) // Match multiple characters (varying length)
			{
				// Check if bound, if not then branch
				if (rule[offset].isBound) 
				{
					if (bound.count(rule[offset].label))
					{
						string tmpStr = bound[rule[offset].label];
						
						tmpStr = L.conjugateString(tmpStr,rule[offset].doConjugate);
						
						for (i=0;(!stop)&&(i<tmpStr.length())&&(i+strpos < matchstr.length());i++)
						{
							if (tmpStr[i] != matchstr[i]) stop = 1; // No match
						}
						
						if (i<tmpStr.length())
						{
							stop = 1; // Didn't finish!
						}
												
						strpos += tmpStr.length();
					}
					else
					{
						branchTree(offset, strpos, rule, matchstr, L);
						for (i=0;i<subTrees.size();i++)
						{
							subTrees[i].applyRule(rule,matchstr,L);
						}
						return;
					}
				}
				else
				{
					branchTree(offset, strpos, rule, matchstr, L);
					for (i=0;i<subTrees.size();i++)
					{
						subTrees[i].applyRule(rule,matchstr,L);
					}
					return;
				}
			}
			else
			{
				if (!isMatch(rule[offset], matchstr[strpos],L))
				{
					stop = 1; // no match
				}
				
				strpos++;
			}
		}			
		
		offset++; 		
	} while (!stop && (offset<rule.size()));
	
	if (stop)
	{
		// No match!		
		isNull = true;
	}
	else // Got to the end!
	{
		if (strpos < matchstr.length())
		{
			isNull = true; // Did not finish matching compound
		}
		else
			isEnd = true;
	}
}

void SearchSubtree::branchTree(int offset, int strpos, vector<Symbol> rule, string matchstr, Library L)
{
	string boundString = "";
	// Everything here is not yet bound
	
	if (offset == rule.size()-1) // This is the terminal wildcard...
	{
		int totalmatch = 1;
		
		for (int i=0;i<matchstr.length()-strpos;i++) // Try all remaining end-characters for the wildcard
		{
			char C = matchstr[i+strpos];
			
			if (!L.isMemberOfSet(C, rule[offset].label, rule[offset].doConjugate)) totalmatch = 0;
		}
		
		if (!totalmatch)
		{
			isNull = true;
			return;
		}
		else
		{
			isEnd = true;
			if (rule[offset].isBound)
			{
				boundString = matchstr.substr(strpos, matchstr.length()-strpos);
				bound[rule[offset].label] = boundString;
			}
			return;
		}
	}
	else
	{
		char C,C2;
		int i;
		
		// Bit of a hack - match zero-length wildcards
		// Should bundle this elsewhere to avoid code duplication
		i=-1;
		C2 = matchstr[strpos];
		if (L.isMemberOfSet(C2, rule[offset+1].label, rule[offset+1].doConjugate))
		{
			SearchSubtree Child;
				
			Child.bound = bound;
			if (rule[offset].isBound)
			{
				boundString = matchstr.substr(strpos, i+1);
				Child.bound[rule[offset].label] = boundString;
			}
			Child.stringpos = strpos+i+1;
			Child.rulepos = offset+1;
			Child.isNull = false;
			Child.isEnd = false;
			subTrees.push_back(Child);
		}
		
		// Now lets match wildcards of non-zero length
		for (i=0;i<matchstr.length()-strpos;i++) // Try all remaining end-characters for the wildcard
		{
			C = matchstr[i+strpos]; 
			if (i+strpos+1 < matchstr.length())
				C2 = matchstr[i+strpos+1];
				
			if (!L.isMemberOfSet(C, rule[offset].label, rule[offset].doConjugate))
			{
				if (!subTrees.size()) isNull = true; // Found no matches, so this is a dead end
				return; // Mismatch ends the sequence
			}
			else
			{
				if ( ( (i+strpos+1 == matchstr.length()) && (rule[offset+1].isWild) && (!bound.count(rule[offset+1].label)) ) 
					 || L.isMemberOfSet(C2, rule[offset+1].label, rule[offset+1].doConjugate))
				{
					SearchSubtree Child;
				
					Child.bound = bound;
					if (rule[offset].isBound)
					{
						boundString = matchstr.substr(strpos, i+1);
						Child.bound[rule[offset].label] = boundString;
					}
					Child.stringpos = strpos+i+1;
					Child.rulepos = offset+1;
					Child.isNull = false;
					Child.isEnd = false;
					subTrees.push_back(Child);
				}
			}
		}		
	}
}

void ReactionRule::parseRule(Library L)
{
	int i, stop = 0;	
	Nreac = 0;
	string substr;
	vector<Symbol> curRule;
	Symbol curSymbol;
	int inrule = 0;
	int prodMode = 0;
	Nprod = 0;

	curSymbol.clear();
	
	for (i=0;(i<rule.length())&&!stop;i++)
	{		
		if (rule[i] == '+') { 
			if (inrule) { curRule.push_back(curSymbol); curSymbol.clear(); inrule = 0; }
			if (!prodMode) { reacRules.push_back(curRule); Nreac++; }
			else { prodRules.push_back(curRule); Nprod++; }
			curRule.clear(); }
		else if (rule[i] == '=') { 
			if (inrule) { curRule.push_back(curSymbol); curSymbol.clear(); inrule = 0; }
			if (!prodMode) { reacRules.push_back(curRule); Nreac++; }
			else { prodRules.push_back(curRule); Nprod++; }
			curRule.clear(); prodMode = 1; }
		else 
		{
			if ((rule[i]!=' ')&&(rule[i]!='\n')) // remove whitespace
			{
				if (!inrule)
				{
					inrule = 1;
					curSymbol.label = rule[i];
					if (find(L.atomTypes.begin(),L.atomTypes.end(),rule[i])!=L.atomTypes.end())
					{
						curSymbol.isSingleton = true;
					}
					else
					{
						curSymbol.isBound = true;
					}
				}
				else
				{
					switch (rule[i])
					{
						case '&':
								i++;
								curSymbol.doConjugate=rule[i]-'0'; 
							break; // Conjugate modifier
						case '\'':
							curSymbol.isReversed = true;
							break; // reverse modifier
						case '*':
							curSymbol.isWild = true;
							break; // wildcard modifier
						default: // Some other character, so go to the next rule and reparse this character
							if (inrule) { curRule.push_back(curSymbol); curSymbol.clear(); inrule = 0; i--;}
							break;
					}
				}
			} 
			else
			{
				if (inrule) { curRule.push_back(curSymbol); curSymbol.clear(); inrule = 0; }
			}
		}
	}
	
	if (inrule) { curRule.push_back(curSymbol); curSymbol.clear(); inrule = 0; }
	if (!prodMode) { reacRules.push_back(curRule); Nreac++; }
	else { prodRules.push_back(curRule); Nprod++; }
}

vector<vector<string>> ReactionRule::getAllProducts(vector<string> reactants, Library L)
{
	SearchSubtree Root;
	int i=0;
	
	Root.isNull = false; Root.isEnd = true;
	Root.rulepos = 0; Root.stringpos = 0;
	
	for (i=0;i<reacRules.size();i++)
		Root.parseOnLeaves(reacRules[i], reactants[i], L);	
	
	return Root.generateProducts(prodRules, L);		
}

string writeRule(vector<Symbol> rule)
{
	int i;
	string rstr = "";
	
	for (i=0;i<rule.size();i++)
	{
		rstr = rstr + rule[i].label;
		if (rule[i].isWild) rstr = rstr + "*";
		if (rule[i].isReversed) rstr = rstr + "'";
		if (rule[i].doConjugate) rstr = rstr + "&" + to_string(rule[i].doConjugate);
	}
	
	return rstr;
}
