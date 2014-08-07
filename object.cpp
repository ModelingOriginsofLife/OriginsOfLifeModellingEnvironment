Object::Object()
{
}

Object::Object(vector<DataType*> schema)
{
	for (int i=0;i<schema.size();i++)
	{
		DataType *D = schema[i]->newMember();
		
		properties.push_back(D);
	}
}
