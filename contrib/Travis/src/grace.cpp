/*****************************************************************************
    TRAVIS - Trajectory Analyzer and Visualizer
    http://www.travis-analyzer.de/

    Copyright (c) 2009-2016 Martin Brehm
                  2012-2016 Martin Thomas

    This file written by Martin Brehm.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include "grace.h"
#include <math.h>
#include <stdlib.h>
#include "tools.h"
#include "maintools.h"


CGrace::CGrace()
{
	m_iBGColor = 0xFFFFFF;
	AddGraph();

	AddColor(255,255,255,"white");
	AddColor(  0,  0,  0,"black");
	AddColor(255,  0,  0,"red");
	AddColor(  0,255,  0,"green");
	AddColor(  0,  0,255,"blue");
	AddColor(255,255,  0,"yellow");
	AddColor(188,143,143,"brown");
	AddColor(220,220,220,"grey");
	AddColor(148,  0,211,"violet");
	AddColor(  0,255,255,"cyan");
	AddColor(255,  0,255,"magenta");
	AddColor(255,165,  0,"orange");
	AddColor(114, 33,188,"indigo");
	AddColor(103,  7, 72,"maroon");
	AddColor( 64,224,208,"turquoise");
	AddColor(  0,139,  0,"green4");

	m_oaGraceGraphs.SetName("CGrace::m_oaGraceGraphs");
	m_oaGraceColors.SetName("CGrace::m_oaGraceColors");
}


CGrace::~CGrace()
{

}


CGraceGraph::CGraceGraph()
{
	try { m_sLabelX = new char[1]; } catch(...) { m_sLabelX = NULL; }
	if (m_sLabelX == NULL) NewException((double)sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	*m_sLabelX = 0;

	try { m_sLabelY = new char[1]; } catch(...) { m_sLabelY = NULL; }
	if (m_sLabelY == NULL) NewException((double)sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	*m_sLabelY = 0;

	try { m_sTitle = new char[1]; } catch(...) { m_sTitle = NULL; }
	if (m_sTitle == NULL) NewException((double)sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	*m_sTitle = 0;

	try { m_sSubTitle = new char[1]; } catch(...) { m_sSubTitle = NULL; }
	if (m_sSubTitle == NULL) NewException((double)sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_fYAxisBarWidth = 1.0;
	*m_sSubTitle = 0;
	m_bShowFrame = true;
	m_bShowXAxis = true;
	m_bTicksBothSidesX = false;
	m_bTicksBothSidesY = false;
	m_bTickLabelsBothSidesX = false;
	m_bTickLabelsBothSidesY = false;
	m_bLegend = false;
	m_bTickInX = true;
	m_bTickInY = true;
	m_bTickLabels = true;
	m_bInvertXAxis = false;
	m_bInvertYAxis = false;
	m_bTicks = true;
	m_fFrameWidth = 2.0;
	m_fViewportX1 = 0.15;
	m_fViewportY1 = 0.15;
	m_fViewportX2 = 1.15;
	m_fViewportY2 = 0.85;
	m_bInvert = false;
	m_iTickMinorX = 1;
	m_iTickMinorY = 1;

	m_oaDatasets.SetName("CGraceGraph::m_oaDatasets");
	m_oaCustomLabelsX.SetName("CGraceGraph::m_oaCustomLabelsX");
	m_oaCustomLabelsY.SetName("CGraceGraph::m_oaCustomLabelsY");
	m_oaLines.SetName("CGraceGraph::m_oaLines");
}


CGraceGraph::~CGraceGraph()
{

}


CGraceDataset::CGraceDataset()
{
	try { m_sName = new char[1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	*m_sName = 0;
	m_iLineColorIndex = 0;
	m_iLineStyle = 0;
	m_fLineWidth = 2.0f;
	m_iSymbColorIndex = 0;
	m_iBegin = 0;
	m_iEnd = 0;
	m_bFill = false;

	m_faValues.SetName("CGraceDataset::m_faValues");
}


CGraceDataset::~CGraceDataset()
{
}


CGraceColor::CGraceColor()
{
	try { m_sName = new char[1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	*m_sName = 0;
}


CGraceColor::~CGraceColor()
{

}


void CGrace::MakeTicks()
{
//	float d1, d2;
	CGraceGraph *g;

	g = CurrentGraph();

/*	d1 = dec(m_fMinValX,2);
	d2 = dec(m_fMaxValX,2);
	m_fMinRangeX = minbound(m_fMinValX,d1);
	m_fMaxRangeX = maxbound(m_fMaxValX,d2);*/
	g->m_fTickMajorX = majorticks(g->m_fMinRangeX,g->m_fMaxRangeX);
	g->m_iTickPrecX = (int)max(0,1-int(ceil(log10(g->m_fTickMajorX))));

/*	d1 = dec(m_fMinValY,2);
	d2 = dec(m_fMaxValY,2);
	m_fMinRangeY = minbound(m_fMinValY,d1);
	m_fMaxRangeY = maxbound(m_fMaxValY,d2);*/
	g->m_fTickMajorY = majorticks(g->m_fMinRangeY,g->m_fMaxRangeY);
	g->m_iTickPrecY = (int)max(0,1-int(ceil(log10(g->m_fTickMajorY))));
}


void CGrace::WriteAgr(const char *s, bool silent)
{
	int z, z2/*, i*/;
	FILE *a;
//	CGraceGraph *g;

	if ((!silent) && (CurrentGraph()->m_oaDatasets.GetSize() > 2))
		mprintf("      Writing file...\n");

	a = OpenFileWrite(s,true);

	mfprintf(a,"# Grace project file\n");
	mfprintf(a,"# Written by TRAVIS, the Trajectory Analyser and Visualizer\n");
	mfprintf(a,"# See http://www.travis-analyzer.de\n");
	mfprintf(a,"@version 50119\n");
	mfprintf(a,"@page size 792, 612\n");
	mfprintf(a,"@page scroll 5%c\n",'%');
	mfprintf(a,"@page inout 5%c\n",'%');
	mfprintf(a,"@link page off\n");
	mfprintf(a,"@map font 0 to \"Times-Roman\", \"Times-Roman\"\n");
	mfprintf(a,"@map font 1 to \"Times-Italic\", \"Times-Italic\"\n");
	mfprintf(a,"@map font 2 to \"Times-Bold\", \"Times-Bold\"\n");
	mfprintf(a,"@map font 3 to \"Times-BoldItalic\", \"Times-BoldItalic\"\n");
	mfprintf(a,"@map font 4 to \"Helvetica\", \"Helvetica\"\n");
	mfprintf(a,"@map font 5 to \"Helvetica-Oblique\", \"Helvetica-Oblique\"\n");
	mfprintf(a,"@map font 6 to \"Helvetica-Bold\", \"Helvetica-Bold\"\n");
	mfprintf(a,"@map font 7 to \"Helvetica-BoldOblique\", \"Helvetica-BoldOblique\"\n");
	mfprintf(a,"@map font 8 to \"Courier\", \"Courier\"\n");
	mfprintf(a,"@map font 9 to \"Courier-Oblique\", \"Courier-Oblique\"\n");
	mfprintf(a,"@map font 10 to \"Courier-Bold\", \"Courier-Bold\"\n");
	mfprintf(a,"@map font 11 to \"Courier-BoldOblique\", \"Courier-BoldOblique\"\n");
	mfprintf(a,"@map font 12 to \"Symbol\", \"Symbol\"\n");
	mfprintf(a,"@map font 13 to \"ZapfDingbats\", \"ZapfDingbats\"\n");

	for (z=0;z<m_oaGraceColors.GetSize();z++)
		mfprintf(a,"@map color %d to (%d, %d, %d), \"%s\"\n",z,((CGraceColor*)m_oaGraceColors[z])->m_iColorRed,((CGraceColor*)m_oaGraceColors[z])->m_iColorGreen,((CGraceColor*)m_oaGraceColors[z])->m_iColorBlue,((CGraceColor*)m_oaGraceColors[z])->m_sName);

	mfprintf(a,"@map color %d to (%lu, %lu, %lu), \"color_bg\"\n",m_oaGraceColors.GetSize(),m_iBGColor % 256,(m_iBGColor >> 8) % 256,(m_iBGColor >> 16) % 256);

/*	mfprintf(a,"@map color 1 to (0, 0, 0), \"black\"\n");
	mfprintf(a,"@map color 2 to (255, 0, 0), \"red\"\n");
	mfprintf(a,"@map color 3 to (0, 255, 0), \"green\"\n");
	mfprintf(a,"@map color 4 to (0, 0, 255), \"blue\"\n");
	mfprintf(a,"@map color 5 to (255, 255, 0), \"yellow\"\n");
	mfprintf(a,"@map color 6 to (188, 143, 143), \"brown\"\n");
	mfprintf(a,"@map color 7 to (220, 220, 220), \"grey\"\n");
	mfprintf(a,"@map color 8 to (148, 0, 211), \"violet\"\n");
	mfprintf(a,"@map color 9 to (0, 255, 255), \"cyan\"\n");
	mfprintf(a,"@map color 10 to (255, 0, 255), \"magenta\"\n");
	mfprintf(a,"@map color 11 to (255, 165, 0), \"orange\"\n");
	mfprintf(a,"@map color 12 to (114, 33, 188), \"indigo\"\n");
	mfprintf(a,"@map color 13 to (103, 7, 72), \"maroon\"\n");
	mfprintf(a,"@map color 14 to (64, 224, 208), \"turquoise\"\n");
	mfprintf(a,"@map color 15 to (0, 139, 0), \"green4\"\n");

	i = 16;
	for (z=0;z<m_oaGraceGraphs.GetSize();z++)
	{
		g = (CGraceGraph*)m_oaGraceGraphs[z];
		g->m_iColorIndex = i;
		for (z2=0;z2<g->m_oaDatasets.GetSize();z2++)
		{
			mfprintf(a,"@map color %d to (%lu, %lu, %lu), \"color_line%d_%d\"\n",z2*2+g->m_iColorIndex,((CGraceDataset*)g->m_oaDatasets[z2])->m_iLineColor % 256,(((CGraceDataset*)g->m_oaDatasets[z2])->m_iLineColor >> 8) % 256,(((CGraceDataset*)g->m_oaDatasets[z2])->m_iLineColor >> 16) % 256,z+1,z2+1);
			mfprintf(a,"@map color %d to (%lu, %lu, %lu), \"color_symb%d_%d\"\n",z2*2+g->m_iColorIndex+1,((CGraceDataset*)g->m_oaDatasets[z2])->m_iSymbColor % 256,(((CGraceDataset*)g->m_oaDatasets[z2])->m_iSymbColor >> 8) % 256,(((CGraceDataset*)g->m_oaDatasets[z2])->m_iSymbColor >> 16) % 256,z+1,z2+1);
			i+=2;
		}
	}

	mfprintf(a,"@map color %d to (%lu, %lu, %lu), \"color_bg\"\n",i,m_iBGColor % 256,(m_iBGColor >> 8) % 256,(m_iBGColor >> 16) % 256);
*/
	mfprintf(a,"@reference date 0\n");
	mfprintf(a,"@date wrap off\n");
	mfprintf(a,"@date wrap year 1950\n");
	mfprintf(a,"@default linewidth 1.0\n");
	mfprintf(a,"@default linestyle 1\n");
	mfprintf(a,"@default color 1\n");
	mfprintf(a,"@default pattern 1\n");
	mfprintf(a,"@default font 0\n");
	mfprintf(a,"@default char size 1.000000\n");
	mfprintf(a,"@default symbol size 1.000000\n");
	mfprintf(a,"@default sformat \"%c.8g\"\n",'%');
	mfprintf(a,"@background color %d\n",m_oaGraceColors.GetSize());
	mfprintf(a,"@page background fill on\n");
	mfprintf(a,"@timestamp off\n");
	mfprintf(a,"@timestamp 0.03, 0.03\n");
	mfprintf(a,"@timestamp color 1\n");
	mfprintf(a,"@timestamp rot 0\n");
	mfprintf(a,"@timestamp font 0\n");
	mfprintf(a,"@timestamp char size 1.000000\n");
	mfprintf(a,"@timestamp def \"xx\"\n");

	for (z=0;z<m_oaGraceGraphs.GetSize();z++)
		for (z2=0;z2<((CGraceGraph*)m_oaGraceGraphs[z])->m_oaLines.GetSize();z2++)
			((CGraceLine*)((CGraceGraph*)m_oaGraceGraphs[z])->m_oaLines[z2])->WriteLine(a,((CGraceGraph*)m_oaGraceGraphs[z])->m_iNumber);

	mfprintf(a,"@r0 off\n");
	mfprintf(a,"@link r0 to g0\n");
	mfprintf(a,"@r0 type above\n");
	mfprintf(a,"@r0 linestyle 1\n");
	mfprintf(a,"@r0 linewidth 1.0\n");
	mfprintf(a,"@r0 color 1\n");
	mfprintf(a,"@r0 line 0, 0, 0, 0\n");
	mfprintf(a,"@r1 off\n");
	mfprintf(a,"@link r1 to g0\n");
	mfprintf(a,"@r1 type above\n");
	mfprintf(a,"@r1 linestyle 1\n");
	mfprintf(a,"@r1 linewidth 1.0\n");
	mfprintf(a,"@r1 color 1\n");
	mfprintf(a,"@r1 line 0, 0, 0, 0\n");
	mfprintf(a,"@r2 off\n");
	mfprintf(a,"@link r2 to g0\n");
	mfprintf(a,"@r2 type above\n");
	mfprintf(a,"@r2 linestyle 1\n");
	mfprintf(a,"@r2 linewidth 1.0\n");
	mfprintf(a,"@r2 color 1\n");
	mfprintf(a,"@r2 line 0, 0, 0, 0\n");
	mfprintf(a,"@r3 off\n");
	mfprintf(a,"@link r3 to g0\n");
	mfprintf(a,"@r3 type above\n");
	mfprintf(a,"@r3 linestyle 1\n");
	mfprintf(a,"@r3 linewidth 1.0\n");
	mfprintf(a,"@r3 color 1\n");
	mfprintf(a,"@r3 line 0, 0, 0, 0\n");
	mfprintf(a,"@r4 off\n");
	mfprintf(a,"@link r4 to g0\n");
	mfprintf(a,"@r4 type above\n");
	mfprintf(a,"@r4 linestyle 1\n");
	mfprintf(a,"@r4 linewidth 1.0\n");
	mfprintf(a,"@r4 color 1\n");
	mfprintf(a,"@r4 line 0, 0, 0, 0\n");

	for (z=0;z<m_oaGraceGraphs.GetSize();z++)
		((CGraceGraph*)m_oaGraceGraphs[z])->Write(a);

	for (z=0;z<m_oaGraceGraphs.GetSize();z++)
		((CGraceGraph*)m_oaGraceGraphs[z])->WriteData(a,silent);

	fclose(a);
}


void CGraceGraph::WriteData(FILE *a, bool silent)
{
	int z;

	if ((!silent) && (m_oaDatasets.GetSize() > 2))
		mprintf(WHITE,"        [");
	for (z=0;z<m_oaDatasets.GetSize();z++)
	{
		mfprintf(a,"@target G%d.S%d\n",m_iNumber,((CGraceDataset*)m_oaDatasets[z])->m_iNumber);
		((CGraceDataset*)m_oaDatasets[z])->WriteSet(a);
		if ((!silent) && (m_oaDatasets.GetSize() > 2))
			if ((z % max(1,(m_oaDatasets.GetSize()/40))) == 0)
				mprintf(WHITE,"#");
	}
	if ((!silent) && (m_oaDatasets.GetSize() > 2))
		mprintf(WHITE,"]\n");
}


CGraceDataset* CGraceGraph::CurrentDataset()
{
	return (CGraceDataset*)m_oaDatasets[m_oaDatasets.GetSize()-1];
}


CGraceDataset* CGraceGraph::Dataset(int i)
{
	if ((i >= m_oaDatasets.GetSize()) || (i < 0))
	{
		eprintf("CGraceGraph::Dataset(): %d/%d.\n",i,m_oaDatasets.GetSize());
		abort();
	}
	return (CGraceDataset*)m_oaDatasets[m_oaDatasets.GetSize()-1];
}


void CGraceGraph::Write(FILE *a)
{
	int z;

	mfprintf(a,"@g%d on\n",m_iNumber);
	mfprintf(a,"@g%d hidden false\n",m_iNumber);
	mfprintf(a,"@g%d type XY\n",m_iNumber);
	mfprintf(a,"@g%d stacked false\n",m_iNumber);
	mfprintf(a,"@g%d bar hgap 0.000000\n",m_iNumber);
	mfprintf(a,"@g%d fixedpoint off\n",m_iNumber);
	mfprintf(a,"@g%d fixedpoint type 0\n",m_iNumber);
	mfprintf(a,"@g%d fixedpoint xy 0.000000, 0.000000\n",m_iNumber);
	mfprintf(a,"@g%d fixedpoint format general general\n",m_iNumber);
	mfprintf(a,"@g%d fixedpoint prec 6, 6\n",m_iNumber);
	mfprintf(a,"@with g%d\n",m_iNumber);
	mfprintf(a,"@    world %f, %f, %f, %f\n",m_fMinRangeX,m_fMinRangeY,m_fMaxRangeX,m_fMaxRangeY);
	mfprintf(a,"@    stack world 0, 0, 0, 0\n");
	mfprintf(a,"@    znorm 1\n");
	mfprintf(a,"@    view %f, %f, %f, %f\n",m_fViewportX1,m_fViewportY1,m_fViewportX2,m_fViewportY2);
	mfprintf(a,"@    title \"%s\"\n",m_sTitle);
	mfprintf(a,"@    title font 0\n");
	mfprintf(a,"@    title size 1.500000\n");
	mfprintf(a,"@    title color 1\n");
	mfprintf(a,"@    subtitle \"%s\"\n",m_sSubTitle);
	mfprintf(a,"@    subtitle font 0\n");
	mfprintf(a,"@    subtitle size 1.000000\n");
	mfprintf(a,"@    subtitle color 1\n");
	mfprintf(a,"@    xaxes scale Normal\n");
	mfprintf(a,"@    yaxes scale Normal\n");
	if (m_bInvertXAxis)
		mfprintf(a,"@    xaxes invert on\n");
			else mfprintf(a,"@    xaxes invert off\n");
	if (m_bInvertYAxis)
		mfprintf(a,"@    yaxes invert on\n");
			else mfprintf(a,"@    yaxes invert off\n");
	mfprintf(a,"@    xaxis  on\n");
	mfprintf(a,"@    xaxis  type zero false\n");
	mfprintf(a,"@    xaxis  offset 0.000000 , 0.000000\n");
	mfprintf(a,"@    xaxis  bar %s\n",m_bShowXAxis?"on":"off");
	mfprintf(a,"@    xaxis  bar color 1\n");
	mfprintf(a,"@    xaxis  bar linestyle 1\n");
	mfprintf(a,"@    xaxis  bar linewidth 1.0\n");
	mfprintf(a,"@    xaxis  label \"%s\"\n",m_bShowXAxis?m_sLabelX:"");
	mfprintf(a,"@    xaxis  label layout para\n");
	mfprintf(a,"@    xaxis  label place auto\n");
	mfprintf(a,"@    xaxis  label char size 1.000000\n");
	mfprintf(a,"@    xaxis  label font 4\n");
	mfprintf(a,"@    xaxis  label color 1\n");
	mfprintf(a,"@    xaxis  label place normal\n");
	if (m_bTicks)
		mfprintf(a,"@    xaxis  tick %s\n",m_bShowXAxis?"on":"off");
			else mfprintf(a,"@    xaxis  tick off\n");
	mfprintf(a,"@    xaxis  tick major %f\n",m_fTickMajorX);
	mfprintf(a,"@    xaxis  tick minor ticks %d\n",m_iTickMinorX);
	mfprintf(a,"@    xaxis  tick default 6\n");
	mfprintf(a,"@    xaxis  tick place rounded true\n");
	if (m_bTickInX)
		mfprintf(a,"@    xaxis  tick in\n");
			else mfprintf(a,"@    xaxis  tick out\n");
	mfprintf(a,"@    xaxis  tick major size 1.000000\n");
	mfprintf(a,"@    xaxis  tick major color 1\n");
	mfprintf(a,"@    xaxis  tick major linewidth 2.0\n");
	mfprintf(a,"@    xaxis  tick major linestyle 1\n");
	mfprintf(a,"@    xaxis  tick major grid off\n");
	mfprintf(a,"@    xaxis  tick minor color 1\n");
	mfprintf(a,"@    xaxis  tick minor linewidth 2.0\n");
	mfprintf(a,"@    xaxis  tick minor linestyle 1\n");
	mfprintf(a,"@    xaxis  tick minor grid off\n");
	mfprintf(a,"@    xaxis  tick minor size 0.500000\n");
	if (m_bTickLabels)
		mfprintf(a,"@    xaxis  ticklabel %s\n",m_bShowXAxis?"on":"off");
			else mfprintf(a,"@    xaxis  ticklabel off\n");
	mfprintf(a,"@    xaxis  ticklabel format decimal\n");
	mfprintf(a,"@    xaxis  ticklabel prec %d\n",m_iTickPrecX);
	mfprintf(a,"@    xaxis  ticklabel formula \"\"\n");
	mfprintf(a,"@    xaxis  ticklabel append \"\"\n");
	mfprintf(a,"@    xaxis  ticklabel prepend \"\"\n");
	mfprintf(a,"@    xaxis  ticklabel angle 0\n");
	mfprintf(a,"@    xaxis  ticklabel skip 0\n");
	mfprintf(a,"@    xaxis  ticklabel stagger 0\n");
	if (m_bTickLabelsBothSidesX)
		mfprintf(a,"@    xaxis  ticklabel place both\n");
			else mfprintf(a,"@    xaxis  ticklabel place normal\n");
	mfprintf(a,"@    xaxis  ticklabel offset auto\n");
	mfprintf(a,"@    xaxis  ticklabel offset 0.000000 , 0.010000\n");
	mfprintf(a,"@    xaxis  ticklabel start type auto\n");
	mfprintf(a,"@    xaxis  ticklabel start 0.000000\n");
	mfprintf(a,"@    xaxis  ticklabel stop type auto\n");
	mfprintf(a,"@    xaxis  ticklabel stop 0.000000\n");
	mfprintf(a,"@    xaxis  ticklabel char size 1.000000\n");
	mfprintf(a,"@    xaxis  ticklabel font 4\n");
	mfprintf(a,"@    xaxis  ticklabel color 1\n");
	if (m_bTicksBothSidesX)
		mfprintf(a,"@    xaxis  tick place both\n");
			else mfprintf(a,"@    xaxis  tick place normal\n");
	if (m_oaCustomLabelsX.GetSize() > 0)
	{
		mfprintf(a,"@    xaxis  tick spec type both\n");
		mfprintf(a,"@    xaxis  tick spec %d\n",m_oaCustomLabelsX.GetSize());
		for (z=0;z<m_oaCustomLabelsX.GetSize();z++)
			((CGraceCustomLabel*)m_oaCustomLabelsX[z])->Write(a);
	} else
	{
		mfprintf(a,"@    xaxis  tick spec type none\n");
	}
	mfprintf(a,"@    yaxis  on\n");
	mfprintf(a,"@    yaxis  type zero false\n");
	mfprintf(a,"@    yaxis  offset 0.000000 , 0.000000\n");
	mfprintf(a,"@    yaxis  bar on\n");
	mfprintf(a,"@    yaxis  bar color 1\n");
	mfprintf(a,"@    yaxis  bar linestyle 1\n");
	mfprintf(a,"@    yaxis  bar linewidth %.1f\n",m_fYAxisBarWidth);
	mfprintf(a,"@    yaxis  label \"%s\"\n",m_sLabelY);
	mfprintf(a,"@    yaxis  label layout para\n");
	mfprintf(a,"@    yaxis  label place auto\n");
	mfprintf(a,"@    yaxis  label char size 1.000000\n");
	mfprintf(a,"@    yaxis  label font 4\n");
	mfprintf(a,"@    yaxis  label color 1\n");
	mfprintf(a,"@    yaxis  label place normal\n");
	if (m_bTicks)
		mfprintf(a,"@    yaxis  tick on\n");
			else mfprintf(a,"@    yaxis  tick off\n");
	mfprintf(a,"@    yaxis  tick major %f\n",m_fTickMajorY);
	mfprintf(a,"@    yaxis  tick minor ticks %d\n",m_iTickMinorY);
	mfprintf(a,"@    yaxis  tick default 6\n");
	mfprintf(a,"@    yaxis  tick place rounded true\n");
	if (m_bTickInY)
		mfprintf(a,"@    yaxis  tick in\n");
			else mfprintf(a,"@    yaxis  tick out\n");
	mfprintf(a,"@    yaxis  tick major size 1.000000\n");
	mfprintf(a,"@    yaxis  tick major color 1\n");
	mfprintf(a,"@    yaxis  tick major linewidth 2.0\n");
	mfprintf(a,"@    yaxis  tick major linestyle 1\n");
	mfprintf(a,"@    yaxis  tick major grid off\n");
	mfprintf(a,"@    yaxis  tick minor color 1\n");
	mfprintf(a,"@    yaxis  tick minor linewidth 2.0\n");
	mfprintf(a,"@    yaxis  tick minor linestyle 1\n");
	mfprintf(a,"@    yaxis  tick minor grid off\n");
	mfprintf(a,"@    yaxis  tick minor size 0.500000\n");
	if (m_bTickLabels)
		mfprintf(a,"@    yaxis  ticklabel on\n");
			else mfprintf(a,"@    yaxis  ticklabel off\n");
	mfprintf(a,"@    yaxis  ticklabel format decimal\n");
	mfprintf(a,"@    yaxis  ticklabel prec %d\n",m_iTickPrecY);
	mfprintf(a,"@    yaxis  ticklabel formula \"\"\n");
	mfprintf(a,"@    yaxis  ticklabel append \"\"\n");
	mfprintf(a,"@    yaxis  ticklabel prepend \"\"\n");
	mfprintf(a,"@    yaxis  ticklabel angle 0\n");
	mfprintf(a,"@    yaxis  ticklabel skip 0\n");
	mfprintf(a,"@    yaxis  ticklabel stagger 0\n");
	if (m_bTickLabelsBothSidesY)
		mfprintf(a,"@    yaxis  ticklabel place both\n");
			else mfprintf(a,"@    yaxis  ticklabel place normal\n");
	mfprintf(a,"@    yaxis  ticklabel offset auto\n");
	mfprintf(a,"@    yaxis  ticklabel offset 0.000000 , 0.010000\n");
	mfprintf(a,"@    yaxis  ticklabel start type auto\n");
	mfprintf(a,"@    yaxis  ticklabel start 0.000000\n");
	mfprintf(a,"@    yaxis  ticklabel stop type auto\n");
	mfprintf(a,"@    yaxis  ticklabel stop 0.000000\n");
	mfprintf(a,"@    yaxis  ticklabel char size 1.000000\n");
	mfprintf(a,"@    yaxis  ticklabel font 4\n");
	mfprintf(a,"@    yaxis  ticklabel color 1\n");
	if (m_bTicksBothSidesY)
		mfprintf(a,"@    yaxis  tick place both\n");
			else mfprintf(a,"@    yaxis  tick place normal\n");
	if (m_oaCustomLabelsY.GetSize() > 0)
	{
		mfprintf(a,"@    yaxis  tick spec type both\n");
		mfprintf(a,"@    yaxis  tick spec %d\n",m_oaCustomLabelsY.GetSize());
		for (z=0;z<m_oaCustomLabelsY.GetSize();z++)
			((CGraceCustomLabel*)m_oaCustomLabelsY[z])->Write(a);
	} else
	{
		mfprintf(a,"@    yaxis  tick spec type none\n");
	}
	mfprintf(a,"@    altxaxis  off\n");
	mfprintf(a,"@    altyaxis  off\n");
	if (m_bLegend)
		mfprintf(a,"@    legend on\n");
			else mfprintf(a,"@    legend off\n");
	mfprintf(a,"@    legend loctype view\n");
	mfprintf(a,"@    legend 0.9, 0.84\n");
	mfprintf(a,"@    legend box color 1\n");
	mfprintf(a,"@    legend box pattern 1\n");
	mfprintf(a,"@    legend box linewidth 2.0\n");
	mfprintf(a,"@    legend box linestyle 1\n");
	mfprintf(a,"@    legend box fill color 0\n");
	mfprintf(a,"@    legend box fill pattern 1\n");
	mfprintf(a,"@    legend font 4\n");
	mfprintf(a,"@    legend char size 1.000000\n");
	mfprintf(a,"@    legend color 1\n");
	mfprintf(a,"@    legend length 4\n");
	mfprintf(a,"@    legend vgap 1\n");
	mfprintf(a,"@    legend hgap 1\n");
	if (m_bInvert)
		mfprintf(a,"@    legend invert true\n");
			else mfprintf(a,"@    legend invert false\n");
	mfprintf(a,"@    frame type 0\n");
	mfprintf(a,"@    frame linestyle 1\n");
	mfprintf(a,"@    frame linewidth %f\n",m_fFrameWidth);
	mfprintf(a,"@    frame color 1\n");
	mfprintf(a,"@    frame pattern %d\n",m_bShowFrame?1:0);
	mfprintf(a,"@    frame background color 0\n");
	mfprintf(a,"@    frame background pattern 0\n");

	for (z=0;z<m_oaDatasets.GetSize();z++)
		((CGraceDataset*)m_oaDatasets[z])->WriteHeader(a);
}


void CGraceDataset::WriteHeader(FILE *a)
{
//	CGraceGraph *g;

//	g = m_pGraph;
	mfprintf(a,"@    s%d hidden false\n",m_iNumber);
	mfprintf(a,"@    s%d type xy\n",m_iNumber);
	mfprintf(a,"@    s%d symbol 0\n",m_iNumber);
	mfprintf(a,"@    s%d symbol size 1.000000\n",m_iNumber);
	mfprintf(a,"@    s%d symbol color %d\n",m_iNumber,m_iSymbColorIndex);
	mfprintf(a,"@    s%d symbol pattern 1\n",m_iNumber);
	mfprintf(a,"@    s%d symbol fill color %d\n",m_iNumber,m_iSymbColorIndex);
	mfprintf(a,"@    s%d symbol fill pattern 0\n",m_iNumber);
	mfprintf(a,"@    s%d symbol linewidth 1.0\n",m_iNumber);
	mfprintf(a,"@    s%d symbol linestyle 1\n",m_iNumber);
	mfprintf(a,"@    s%d symbol char 65\n",m_iNumber);
	mfprintf(a,"@    s%d symbol char font 0\n",m_iNumber);
	mfprintf(a,"@    s%d symbol skip 0\n",m_iNumber);
	mfprintf(a,"@    s%d line type 1\n",m_iNumber);
	mfprintf(a,"@    s%d line linestyle 1\n",m_iNumber);
	mfprintf(a,"@    s%d line linewidth %.1f\n",m_iNumber,m_fLineWidth);
	mfprintf(a,"@    s%d line color %d\n",m_iNumber,m_iLineColorIndex);
	mfprintf(a,"@    s%d line pattern 1\n",m_iNumber);
	mfprintf(a,"@    s%d baseline type 0\n",m_iNumber);
	mfprintf(a,"@    s%d baseline off\n",m_iNumber);
	mfprintf(a,"@    s%d dropline off\n",m_iNumber);
	if (m_bFill)
		mfprintf(a,"@    s%d fill type 2\n",m_iNumber);
			else mfprintf(a,"@    s%d fill type 0\n",m_iNumber);
	mfprintf(a,"@    s%d fill rule 0\n",m_iNumber);
	mfprintf(a,"@    s%d fill color %d\n",m_iNumber,m_iLineColorIndex);
	mfprintf(a,"@    s%d fill pattern 1\n",m_iNumber);
	mfprintf(a,"@    s%d avalue off\n",m_iNumber);
	mfprintf(a,"@    s%d avalue type 2\n",m_iNumber);
	mfprintf(a,"@    s%d avalue char size 1.000000\n",m_iNumber);
	mfprintf(a,"@    s%d avalue font 0\n",m_iNumber);
	mfprintf(a,"@    s%d avalue color 1\n",m_iNumber);
	mfprintf(a,"@    s%d avalue rot 0\n",m_iNumber);
	mfprintf(a,"@    s%d avalue format general\n",m_iNumber);
	mfprintf(a,"@    s%d avalue prec 3\n",m_iNumber);
	mfprintf(a,"@    s%d avalue prepend \"\"\n",m_iNumber);
	mfprintf(a,"@    s%d avalue append \"\"\n",m_iNumber);
	mfprintf(a,"@    s%d avalue offset 0.000000 , 0.000000\n",m_iNumber);
	mfprintf(a,"@    s%d errorbar on\n",m_iNumber);
	mfprintf(a,"@    s%d errorbar place both\n",m_iNumber);
	mfprintf(a,"@    s%d errorbar color 4\n",m_iNumber);
	mfprintf(a,"@    s%d errorbar pattern 1\n",m_iNumber);
	mfprintf(a,"@    s%d errorbar size 1.000000\n",m_iNumber);
	mfprintf(a,"@    s%d errorbar linewidth 1.0\n",m_iNumber);
	mfprintf(a,"@    s%d errorbar linestyle 1\n",m_iNumber);
	mfprintf(a,"@    s%d errorbar riser linewidth 1.0\n",m_iNumber);
	mfprintf(a,"@    s%d errorbar riser linestyle 1\n",m_iNumber);
	mfprintf(a,"@    s%d errorbar riser clip off\n",m_iNumber);
	mfprintf(a,"@    s%d errorbar riser clip length 0.100000\n",m_iNumber);
	mfprintf(a,"@    s%d comment \"%s\"\n",m_iNumber,m_sName);
	mfprintf(a,"@    s%d legend  \"%s\"\n",m_iNumber,m_sName);
}


void CGraceDataset::WriteSet(FILE *a)
{
	int z;

	mfprintf(a,"@type xy\n");

	for (z=m_iBegin;z<(((m_faValues.GetSize()/2)<m_iEnd)?(m_faValues.GetSize()/2):m_iEnd);z++)
		mfprintf(a,"%G %G\n",m_faValues[z*2],m_faValues[z*2+1]);
	mfprintf(a,"&\n");
}


void CGrace::AddDataset()
{
	CGraceDataset *gd;
	CGraceGraph *g;

	g = CurrentGraph();

	try { gd = new CGraceDataset(); } catch(...) { gd = NULL; }
	if (gd == NULL) NewException((double)sizeof(CGraceDataset),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	gd->m_pGraph = g;
	gd->m_iNumber = g->m_oaDatasets.GetSize();
	gd->m_iLineColorIndex = (gd->m_iNumber % 15) + 1; // alles ausser weiss
	g->m_oaDatasets.Add(gd);
}


void CGrace::AddGraph()
{
	CGraceGraph *g;

	try { g = new CGraceGraph(); } catch(...) { g = NULL; }
	if (g == NULL) NewException((double)sizeof(CGraceGraph),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	g->m_iNumber = m_oaGraceGraphs.GetSize();
	m_oaGraceGraphs.Add(g);
}


void CGrace::AddXYTupel(int set, double x, double y)
{
	CGraceGraph *g;

	g = CurrentGraph();
/*	if (set >= g->m_oaDatasets.GetSize())
	{
		eprintf("CGrace::AddXYTupel(): Set %d >= %d.\n",set,g->m_oaDatasets.GetSize());
		return;
	}*/
	((CGraceDataset*)g->m_oaDatasets[set])->m_faValues.Add(x);
	((CGraceDataset*)g->m_oaDatasets[set])->m_faValues.Add(y);
	((CGraceDataset*)g->m_oaDatasets[set])->m_iEnd = ((CGraceDataset*)g->m_oaDatasets[set])->m_faValues.GetSize()/2;
}


void CGrace::AddXYTupel(double x, double y)
{
	AddXYTupel(CurrentGraph()->m_oaDatasets.GetSize()-1,x,y);
}


void CGrace::FindMinMaxVal()
{
	CGraceDataset *gd;
	int z, z2;
	CGraceGraph *g;

	g = CurrentGraph();

	g->m_fMinValX = (float)9E20;
	g->m_fMinValY = (float)9E20;
	g->m_fMaxValX = (float)-9E20;
	g->m_fMaxValY = (float)-9E20;

	for (z=0;z<g->m_oaDatasets.GetSize();z++)
	{
		gd = (CGraceDataset*)g->m_oaDatasets[z];
		for (z2=0;z2<gd->m_faValues.GetSize()/2;z2++)
		{
			if (gd->m_faValues[z2*2] > g->m_fMaxValX)
				g->m_fMaxValX = (float)gd->m_faValues[z2*2];
			if (gd->m_faValues[z2*2] < g->m_fMinValX)
				g->m_fMinValX = (float)gd->m_faValues[z2*2];
			if (gd->m_faValues[z2*2+1] > g->m_fMaxValY)
				g->m_fMaxValY = (float)gd->m_faValues[z2*2+1];
			if (gd->m_faValues[z2*2+1] < g->m_fMinValY)
				g->m_fMinValY = (float)gd->m_faValues[z2*2+1];
		}
//		mprintf("Dataset %d: X %f - %f, Y %f - %f.\n",z+1,g->m_fMinValX,g->m_fMaxValX,g->m_fMinValY,g->m_fMaxValY);
	}
}


void CGrace::SetTitle(const char *s)
{
	CGraceGraph *g;

	g = CurrentGraph();

	if (g->m_sTitle != NULL)
		delete[] g->m_sTitle;

	try { g->m_sTitle = new char[strlen(s)+1]; } catch(...) { g->m_sTitle = NULL; }
	if (g->m_sTitle == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(g->m_sTitle,s);
}


void CGrace::SetSubTitle(const char *s)
{
	CGraceGraph *g;

	g = CurrentGraph();

	if (g->m_sSubTitle != NULL)
		delete[] g->m_sSubTitle;

	try { g->m_sSubTitle = new char[strlen(s)+1]; } catch(...) { g->m_sSubTitle = NULL; }
	if (g->m_sSubTitle == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(g->m_sSubTitle,s);
}


void CGrace::SetLabelX(const char *s)
{
	CGraceGraph *g;

	g = CurrentGraph();

	if (g->m_sLabelX != NULL)
		delete[] g->m_sLabelX;

	try { g->m_sLabelX = new char[strlen(s)+1]; } catch(...) { g->m_sLabelX = NULL; }
	if (g->m_sLabelX == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(g->m_sLabelX,s);
}


void CGrace::SetLabelY(const char *s)
{
	CGraceGraph *g;

	g = CurrentGraph();

	if (g->m_sLabelY != NULL)
		delete[] g->m_sLabelY;

	try { g->m_sLabelY = new char[strlen(s)+1]; } catch(...) { g->m_sLabelY = NULL; }
	if (g->m_sLabelY == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(g->m_sLabelY,s);
}


void CGrace::SetRangeX(double mi, double ma)
{
	CGraceGraph *g;

	g = CurrentGraph();

	g->m_fMinRangeX = (float)mi;
	g->m_fMaxRangeX = (float)ma;
	g->m_fTickMajorX = majorticks(g->m_fMinRangeX,g->m_fMaxRangeX);
	g->m_iTickPrecX = (int)max(0,1-int(ceil(log10(g->m_fTickMajorX))));
}


void CGrace::SetRangeY(double mi, double ma)
{
	CGraceGraph *g;

	g = CurrentGraph();

	g->m_fMinRangeY = (float)mi;
	g->m_fMaxRangeY = (float)ma;
	g->m_fTickMajorY = majorticks(g->m_fMinRangeY,g->m_fMaxRangeY);
	g->m_iTickPrecY = (int)max(0,1-int(ceil(log10(g->m_fTickMajorY))));
}


void CGrace::DuplicateSet(int set)
{
	CGraceDataset* gd;
	CGraceGraph *g;

	g = CurrentGraph();

	try { gd = new CGraceDataset(); } catch(...) { gd = NULL; }
	if (gd == NULL) NewException((double)sizeof(CGraceDataset),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	gd->CopyFrom((CGraceDataset*)g->m_oaDatasets[set]);
	gd->m_iNumber = g->m_oaDatasets.GetSize();
	g->m_oaDatasets.Add(gd);
}


void CGraceDataset::CopyFrom(CGraceDataset *d)
{
	m_faValues.CopyFrom(&d->m_faValues);
	m_iLineColorIndex = d->m_iLineColorIndex;
	m_fLineWidth = d->m_fLineWidth;
	m_bFill = d->m_bFill;
	m_iSymbColorIndex = d->m_iSymbColorIndex;
	m_pGraph = d->m_pGraph;
}


void CGrace::SetSetRange(int set, int start, int end)
{
	CGraceDataset* gd;
	CGraceGraph *g;

	g = CurrentGraph();

	gd = (CGraceDataset*)g->m_oaDatasets[set];
	gd->m_iBegin = start;
	if (end == -1)
		gd->m_iEnd = gd->m_faValues.GetSize()/2;
			else gd->m_iEnd = end;
}


void CGrace::SetSetLineColor(int set, unsigned char r, unsigned char g, unsigned char b)
{
	int z;
	CGraceColor *gc;
//	char buf[64];
	CxString buf;

	for (z=0;z<m_oaGraceColors.GetSize();z++)
	{
		gc = (CGraceColor*)m_oaGraceColors[z];
		if ((gc->m_iColorRed == r) && (gc->m_iColorGreen == g) && (gc->m_iColorBlue == b))
			goto _done;
	}
	if (m_oaGraceColors.GetSize() < 254)
	{
//		sprintf(buf,"set%d",set+1);
		buf.sprintf("set%d",set+1);
		z = AddColor(r,g,b,buf);
	} else z = 1;
_done:
	((CGraceDataset*)CurrentGraph()->m_oaDatasets[set])->m_iLineColorIndex = z;
	((CGraceDataset*)CurrentGraph()->m_oaDatasets[set])->m_iSymbColorIndex = z;
}


void CGrace::SetSetLineColor(unsigned char r, unsigned char g, unsigned char b)
{
	SetSetLineColor(CurrentGraph()->m_oaDatasets.GetSize()-1,r,g,b);
}


void CGrace::SetSetLineColorLong(int set, unsigned long col)
{
	SetSetLineColor(set,col & 0x100,(col >> 8) & 0x100,col >> 16);
}


void CGrace::SetSetLineWidth(int set, float width)
{
	CGraceDataset* gd;
	CGraceGraph *g;

	g = CurrentGraph();

	gd = (CGraceDataset*)g->m_oaDatasets[set];
	gd->m_fLineWidth = width;
}


void CGraceLine::WriteLine(FILE *a, int graph)
{
	mfprintf(a,"@with line\n");
	mfprintf(a,"@    line on\n");
	mfprintf(a,"@    line loctype world\n");
	mfprintf(a,"@    line g%d\n",graph);
	mfprintf(a,"@    line %f, %f, %f, %f\n",m_fX1,m_fY1,m_fX2,m_fY2);
	mfprintf(a,"@    line linewidth %.1f\n",m_fLineWidth);
	mfprintf(a,"@    line linestyle %d\n",m_iLineStyle);
	mfprintf(a,"@    line color %d\n",m_iLineColorIndex);
	mfprintf(a,"@    line arrow 0\n");
	mfprintf(a,"@    line arrow type 0\n");
	mfprintf(a,"@    line arrow length 1.000000\n");
	mfprintf(a,"@    line arrow layout 1.000000, 1.000000\n");
	mfprintf(a,"@line def\n");
}


void CGrace::AddCustomLabelX(bool major, double val, const char *s)
{
	CGraceCustomLabel *cl;
	CGraceGraph *g;

	g = CurrentGraph();

	try { cl = new CGraceCustomLabel(); } catch(...) { cl = NULL; }
	if (cl == NULL) NewException((double)sizeof(CGraceCustomLabel),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	cl->m_bMajor = major;
	cl->m_fValue = val;

	try { cl->m_sText = new char[strlen(s)+1]; } catch(...) { cl->m_sText = NULL; }
	if (cl->m_sText == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	cl->m_bX = true;
	strcpy(cl->m_sText,s);
	cl->m_iNumber = g->m_oaCustomLabelsX.GetSize();
	g->m_oaCustomLabelsX.Add(cl);
}


void CGrace::AddCustomLabelY(bool major, double val, const char *s)
{
	CGraceCustomLabel *cl;
	CGraceGraph *g;

	g = CurrentGraph();

	try { cl = new CGraceCustomLabel(); } catch(...) { cl = NULL; }
	if (cl == NULL) NewException((double)sizeof(CGraceCustomLabel),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	cl->m_bMajor = major;
	cl->m_fValue = val;

	try { cl->m_sText = new char[strlen(s)+1]; } catch(...) { cl->m_sText = NULL; }
	if (cl->m_sText == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	cl->m_bX = false;
	strcpy(cl->m_sText,s);
	cl->m_iNumber = g->m_oaCustomLabelsY.GetSize();
	g->m_oaCustomLabelsY.Add(cl);
}


void CGrace::AddLine(double x1, double y1, double x2, double y2, double width, int style)
{
	CGraceLine *gl;
	CGraceGraph *g;

	g = CurrentGraph();

	try { gl = new CGraceLine(); } catch(...) { gl = NULL; }
	if (gl == NULL) NewException((double)sizeof(CGraceLine),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	gl->m_fX1 = x1;
	gl->m_fY1 = y1;
	gl->m_fX2 = x2;
	gl->m_fY2 = y2;
	gl->m_fLineWidth = width;
	gl->m_iLineStyle = style;
	gl->m_iLineColorIndex = 1;
	g->m_oaLines.Add(gl);
}


void CGrace::AddLine(double x1, double y1, double x2, double y2, double width, int style, unsigned char r, unsigned char g, unsigned char b)
{
	CGraceLine *gl;
	CGraceGraph *gg;
	CGraceColor *gc;
	int z;
//	char buf[64];
	CxString buf;

	gg = CurrentGraph();

	try { gl = new CGraceLine(); } catch(...) { gl = NULL; }
	if (gl == NULL) NewException((double)sizeof(CGraceLine),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	gl->m_fX1 = x1;
	gl->m_fY1 = y1;
	gl->m_fX2 = x2;
	gl->m_fY2 = y2;
	gl->m_fLineWidth = width;
	gl->m_iLineStyle = style;

	for (z=0;z<m_oaGraceColors.GetSize();z++)
	{
		gc = (CGraceColor*)m_oaGraceColors[z];
		if ((gc->m_iColorRed == r) && (gc->m_iColorGreen == g) && (gc->m_iColorBlue == b))
			goto _done;
	}
	if (m_oaGraceColors.GetSize() < 254)
	{
//		sprintf(buf,"line%d",gg->m_oaLines.GetSize()+1);
		buf.sprintf("line%d",gg->m_oaLines.GetSize()+1);
		z = AddColor(r,g,b,buf);
	} else z = 1;
_done:

	gl->m_iLineColorIndex = z;
	gg->m_oaLines.Add(gl);
}


void CGraceCustomLabel::Write(FILE *a)
{
	if (m_bMajor)
		mfprintf(a,"@    %caxis  tick major %d, %f\n",m_bX?'x':'y',m_iNumber,m_fValue);
			else mfprintf(a,"@    %caxis  tick minor %d, %f\n",m_bX?'x':'y',m_iNumber,m_fValue);
	if (strlen(m_sText) != 0) 
		mfprintf(a,"@    %caxis  ticklabel %d, \"%s\"\n",m_bX?'x':'y',m_iNumber,m_sText);
}


void CGrace::SetViewport(float x1, float y1, float x2, float y2)
{
	CGraceGraph *g;

	g = CurrentGraph();
	g->m_fViewportX1 = x1;
	g->m_fViewportY1 = y1;
	g->m_fViewportX2 = x2;
	g->m_fViewportY2 = y2;
}


CGraceGraph* CGrace::CurrentGraph()
{
	return (CGraceGraph*)m_oaGraceGraphs[m_oaGraceGraphs.GetSize()-1];
}


CGraceDataset* CGrace::LastDataset()
{
	return (CGraceDataset*)CurrentGraph()->m_oaDatasets[CurrentGraph()->m_oaDatasets.GetSize()-1];
}


void CGrace::SetDatasetName(int set, const char *s)
{
	CGraceGraph *g;
	CGraceDataset *ds;

	g = CurrentGraph();

	ds = (CGraceDataset*)g->m_oaDatasets[set];
	if (ds->m_sName != NULL)
		delete[] ds->m_sName;

	try { ds->m_sName = new char[strlen(s)+1]; } catch(...) { ds->m_sName = NULL; }
	if (ds->m_sName == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(ds->m_sName,s);
}


void CGrace::SetDatasetName(const char *s)
{
	SetDatasetName(CurrentGraph()->m_oaDatasets.GetSize()-1,s);
}


void CGrace::WriteCSV(const char *s)
{
	int z, z2;
	FILE *a;
	CGraceGraph *g;

	g = CurrentGraph();

	a = OpenFileWrite(s,true);
	for (z=0;z<((CGraceDataset*)g->m_oaDatasets[0])->m_faValues.GetSize()/2;z++)
	{
		mfprintf(a,"%G;  ",((CGraceDataset*)g->m_oaDatasets[0])->m_faValues[z*2]);
		for (z2=0;z2<g->m_oaDatasets.GetSize();z2++)
		{
			mfprintf(a,"%G",((CGraceDataset*)g->m_oaDatasets[z2])->m_faValues[z*2+1]);
			if (z2 < g->m_oaDatasets.GetSize()-1)
				mfprintf(a,";  ");
		}
		mfprintf(a,"\n");
	}
	fclose(a);
}


int CGrace::AddColor(unsigned char r, unsigned char g, unsigned char b, const char *name)
{
	CGraceColor *gc;

	try { gc = new CGraceColor(); } catch(...) { gc = NULL; }
	if (gc == NULL) NewException((double)sizeof(CGraceColor),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	gc->m_iColorRed = r;
	gc->m_iColorGreen = g;
	gc->m_iColorBlue = b;
	gc->SetName(name);
	m_oaGraceColors.Add(gc);

	return m_oaGraceColors.GetSize()-1;
}


void CGraceColor::SetName(const char *s)
{
	if (m_sName != NULL)
		delete[] m_sName;

	try { m_sName = new char[strlen(s)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sName,s);
}
