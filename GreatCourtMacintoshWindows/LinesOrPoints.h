//Written by Chris J K Williams, University of Bath, UK

void plotLinesOrPoints(int which,int LastVertex)
{
	if(which == 1)glBegin(GL_LINE_STRIP);else glBegin(GL_POINTS);
	int ThisLastVertex = LastVertex;
	if(ThisLastVertex > LastPossibleVertex)ThisLastVertex = LastPossibleVertex;
	for(int myVertex = 0;myVertex <= ThisLastVertex;myVertex ++)
	{
		glColor4f(chris_colours[4 * myVertex],
				  chris_colours[4 * myVertex + 1],
				  chris_colours[4 * myVertex + 2],
				  chris_colours[4 * myVertex + 3]);
		glVertex3f(chris_vertices[3 * myVertex],
				   chris_vertices[3 * myVertex + 1],
				   chris_vertices[3 * myVertex + 2]);
	}
	glEnd();
	//glDrawArrays(GL_LINE_STRIP, 0, LastVertex + 1);
}
