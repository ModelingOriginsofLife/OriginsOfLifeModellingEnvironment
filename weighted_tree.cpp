#include "includes.h"

Node::Node()
{
	isleaf = true;
	child[0] = NULL;
	child[1] = NULL;
	parent = NULL;
}

void Tree::deleteLeaf(Node *N)
{
	if (N==NULL) return;
	if (!N->isleaf) 
	{
		Node *C1 = N->child[0], *C2 = N->child[1]; 
		deleteLeaf(C1); deleteLeaf(C2); return;
	} 

	accessHash.erase(N->value);
	
	Node *parent = N->parent;	
	if (parent == NULL) { delete N; return; }
	
	N->modifyWeight(0.0, 0.0); 
		
	Node *sibling = parent->getOtherChild(N);
	
	// We need to set our grandparent's appropriate child to the sibling rather than moving the sibling's data into this node
	
	if (parent->parent == NULL) // If we have no grandparent, the sibling becomes the root node
	{
		Root = sibling;
		sibling->parent = NULL;
	}
	else
	{
		Node *grandparent = parent->parent;
		if (grandparent->child[0] == parent)
			grandparent->child[0] = sibling;
		else
			grandparent->child[1] = sibling;
			
		sibling->parent = grandparent;			
	}
	
	/*
	 * Parent and one leaf were removed 
	 */
	 
	delete parent;
	delete N;
}

Node *Node::getOtherChild(Node *T)
{
	if (child[0] == T) return child[1];
	else return child[0];
}

Node *Node::findRandomLeaf(WeightingType wType)
{
	if (isleaf) return this;

	/* Bit of an optimization - we only ever use this to find the light branch of the tree, so no need to generate random numbers */
	if (wType == WEIGHT_LIGHT)
	{
		if (child[0]->weight < child[1]->weight) return child[0]->findRandomLeaf(wType);
		else return child[1]->findRandomLeaf(wType);
	}
	
	double r, total;
	
	if (wType == WEIGHT_HEAVY)
	{
		total = child[0]->weight + child[1]->weight;
		r = frand(total);

		if (r<child[0]->weight) return child[0]->findRandomLeaf(wType);
		else return child[1]->findRandomLeaf(wType);
	}
	else if (wType == WEIGHT_HEAVYLENGTH)
	{
		total = child[0]->lweight + child[1]->lweight;
		r = frand(total);

		if (r<child[0]->lweight) return child[0]->findRandomLeaf(wType);
		else return child[1]->findRandomLeaf(wType);
	}	
}

Node *Tree::addNewLeaf(string &val, double w)
{
	if (Root == NULL)
	{
		Root = new Node;
		Root->isleaf = true;
		Root->value = val;
		Root->weight = w;
		Root->lweight = w * val.length();
		Root->parent = NULL;
		Root->child[0] = NULL;
		Root->child[1] = NULL;
		accessHash[val] = Root;
		return Root;
	}
	
	Node *leaf = Root->findRandomLeaf(WEIGHT_LIGHT);
	
	leaf->isleaf = false;
	
	leaf->child[0] = new Node;
	leaf->child[0]->value = leaf->value;
	leaf->child[0]->weight = leaf->weight;
	leaf->child[0]->lweight = leaf->lweight;
	leaf->child[0]->isleaf = true;
	leaf->child[0]->parent = leaf;
	
	accessHash[leaf->child[0]->value] = leaf->child[0];

	leaf->child[1] = new Node;
	leaf->child[1]->value = val;
	leaf->child[1]->weight = w;
	leaf->child[1]->lweight = w * val.length();
	leaf->child[1]->isleaf = true;
	leaf->child[1]->parent = leaf;
	
	accessHash[val] = leaf->child[1];
	
	leaf->modifyWeight(leaf->child[0]->weight + leaf->child[1]->weight, 
					   leaf->child[0]->lweight + leaf->child[1]->lweight);
	
	return leaf->child[1];
}

void Node::modifyWeight(double nw, double nwl)
{
	double dw = nw - weight, 
		   dwl = nwl - lweight;
	
	weight = nw; 
	lweight = nwl;
	
	if (parent!=NULL)
	{
		parent->modifyWeight(parent->weight + dw, 
							 parent->lweight + dwl);
	}
}

void Node::writeContents()
{
	if (isleaf) printf("%.2g %s,",weight,value.c_str());
	else
	{
		child[0]->writeContents();
		child[1]->writeContents();
	}
}

Tree::Tree()
{
	Root = NULL;
}

Tree::~Tree()
{
	if (Root!=NULL)
		deleteLeaf(Root);
	accessHash.clear();
}

void Tree::addSubtree(Node *N)
{
	if (N->isleaf)
	{
		addNewLeaf(N->value, N->weight);
	}
	else
	{
		addSubtree(N->child[0]);
		addSubtree(N->child[1]);
	}
}

Tree &Tree::operator=(const Tree &other)
{
	if (Root != NULL)
		deleteLeaf(Root);
	Root = NULL;
	accessHash.clear();
		
	if (other.Root != NULL)
		addSubtree(other.Root);
	
	return *this;
}
	
Tree::Tree(const Tree &other)
{		
	Root = NULL;
	if (other.Root != NULL)
		addSubtree(other.Root);
}
	
