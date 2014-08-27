#include "includes.h"

void Tree::deleteLeaf(Node *N)
{
	if (N==NULL) return;
	if (!N->isleaf) { deleteLeaf(N->child[0]); deleteLeaf(N->child[1]); deleteLeaf(N); return; } 
	if (N->parent == NULL) { delete N; return; }
	
	N->modifyWeight(0.0); 	
	N->parent->isleaf = true;
	
	accessHash.erase(N->value);
	
	Node *sibling = N->parent->getOtherChild(N);
	
	N->parent->child[0]=N->parent->child[1]=NULL;
	
	N->parent->value = sibling->value;
	N->parent->weight = sibling->weight;
	
	accessHash[N->parent->value] = N->parent;
	
	delete N;
}

Node *Node::getOtherChild(Node *T)
{
	if (child[0] == T) return child[1];
	else return child[0];
}

Node *Node::findRandomLeaf(bool invert)
{
	if (isleaf) return this;

	/* Bit of an optimization - we only ever use this to find the light branch of the tree, so no need to generate random numbers */
	if (invert)
	{
		if (child[0]->weight < child[1]->weight) return child[0]->findRandomLeaf(invert);
		else return child[1]->findRandomLeaf(invert);
	}

	double total = child[0]->weight + child[1]->weight;
	double r = frand(total);
	
	if (r<child[0]->weight) return invert ? child[1]->findRandomLeaf(invert) : child[0]->findRandomLeaf(invert);
	else return invert ? child[0]->findRandomLeaf(invert) : child[1]->findRandomLeaf(invert);
}

Node *Tree::addNewLeaf(string val, double w)
{
	if (Root == NULL)
	{
		Root = new Node;
		Root->isleaf = true;
		Root->value = val;
		Root->weight = w;
		Root->parent = NULL;
		accessHash[val] = Root;
		return Root;
	}
	
	Node *leaf = Root->findRandomLeaf(true);
	
	leaf->isleaf = false;
	leaf->child[0] = new Node;
	leaf->child[0]->value = leaf->value;
	leaf->child[0]->weight = leaf->weight;
	leaf->child[0]->isleaf = true;
	leaf->child[0]->parent = leaf;
	
	accessHash[leaf->child[0]->value] = leaf->child[0];
	
	leaf->child[1] = new Node;
	leaf->child[1]->value = val;
	leaf->child[1]->weight = w;
	leaf->child[1]->isleaf = true;
	leaf->child[1]->parent = leaf;
	
	accessHash[val] = leaf->child[1];
	
	leaf->modifyWeight(leaf->child[0]->weight + leaf->child[1]->weight);
	
	return leaf->child[1];
}

void Node::modifyWeight(double nw)
{
	double dw = nw-weight;
	
	weight = nw;
	
	if (parent!=NULL)
	{
		parent->modifyWeight(parent->weight+dw);
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
