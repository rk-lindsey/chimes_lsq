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

#include "asciiart.h"


CAsciiArt::CAsciiArt()
{
	BTIN;
	m_iPosX = -1;
	m_iPosY = -1;
	m_iSelectedPiece = -1;
	m_oaAsciiPieces.SetName("CAsciiArt::m_oaAsciiPieces");


	AddPiece("                           , \n\
     ,::.::::.            / `. \n\
   ;::        `.         (  ,--. \n\
   ,;           :   _,--'\"\\/\"\"-.\\_ \n\
  ,:             `./ ,---./( O )) `; \n\
  ;  `.           _,'    (  `-' ) /_. \n\
  ;   :         ,'        \\    , (o\\\\ \n\
   ;  :         \\  \\-.  -.__--'   \\' ) \n\
   ;  ;         /\\ (    `-._`-._   \\/ \n\
    ';,        ; : |      -.`._'\\   `. \n\
      ;       ;  : `-.,-..  )  \\'\\   ^. \n\
       ;     ;   `.__   )))\\ ) (`.\\    \\ \n\
        ;   ;        `-`///(, \\ \\ \\)  ,ooo. \n\
         ;,;      ;     ``  ))))_;'(  88888p \n\
          ;      ;         ((-='--',-,Y8888' \n\
          ;     :         ;     ,:'-'  `\"' \n\
           ;     `        ;      | \n\
            ;      (_   __   _,-' \n\
             `---.   ;\"(,-' /                              ____ \n\
      -hrr-       \\ (__/\\\\_`-.-.                     _____/ \n\
                ,(( '/\\/\\/\\`-;;))             ______// \n\
               ((\\''/\\/\\/\\/\\/\\/`/\\      _____/  ____/ \n\
               /'/\\/\\/\\/\\/\\/\\/\\/\\/)  __/  _____/ \n\
              (\\/\\/\\/\\/\\/\\/\\/\\/\\_/ _/ ___/ \n\
               `-|\"\"--\"-.___,--'|-'__/ \n\
                 |              | / \n\
                 |         __,--' \n\
                 _\\,----\"\"'");

	AddPiece("                       ,--\"`-.-._ \n\
                     ,' ,.:::::. `. \n\
                     | ;'     `:: | \n\
                 ,--,-, .::::..;',' \n\
               ,'   ,'      \"\" ,'_ \n\
                '-)\",'o         \"_\\ \n\
                _,' \"\"'  .::.  .:o; \n\
               ;      :::'       ; \n\
               ;     `::;        ; \n\
               : (_;_ ;::        ;\\ \n\
              ;  /   ;::.        ; \\ \n\
             /  /    ,:'        ;:: \\ \n\
          __|__/_   /::        ;   \\ ` \n\
      _,-'       `-<::'___    ;|   \\  \\ \n\
   ,-'    _,-\"\"\";-._`.(__ `.  ;|    \\ |      ____ \n\
  (   ,--'__,--'   |`-`(@)  \\(  `.   `.   ,-'    `-. \n\
   \\___.-'   \\     |::. \\    :    `.   `./,-'\"\"`.   \\ \n\
              \\    |::.  )   : .-.  `-._ ' `--.--'   ) \n\
               \\ .-`.:' /    :      /   `-.__   __,-' \n\
    -hrr-       )    `.'     ;     /         `\"' \n\
               (  `'  ,\\    , ---.( \n\
               ,' --- `:`--'  :  : \\ \n\
              (  :  :  ;   `--`--`-' \n\
               `-`--`--' ");

	AddPiece("                                 ,--.\"\"\n\
                         __,----( o )) \n\
                       ,'--.      , ( \n\
                -\"\",:-(    o ),-'/  ; \n\
                  ( o) `o  _,'\\ / ;( \n\
                   `-;_-<'\\_|-'/ '  ) \n\
                       `.`-.__/ '   | \n\
          \\`.            `. .__,   ; \n\
           )_;--.         \\`       | \n\
          /'(__,-:         )      ; \n\
        ;'    (_,-:     _,::     .| \n\
       ;       ( , ) _,':::'    ,; \n\
      ;         )-,;'  `:'     .:: \n\
      |         `'  ;         `:::\\ \n\
      :       ,'    '            `:\\ \n\
      ;:    '  _,-':         .'     `-. \n\
       ';::..,'  ' ,        `   ,__    `. \n\
         `;''   / ;           _;_,-'     `. \n\
               /            _;--.          \\ \n\
             ,'            / ,'  `.         \\ \n\
            /:            (_(   ,' \\         ) \n\
           /:.               \\_(  /-. .:::,;/ \n\
          (::..                 `-'\\ \"`\"\"' \n\
          ;::::.                    \\        __ \n\
          ,::::::.            .:'    )    ,-'  ) \n\
         /  `;:::::::'`__,:.:::'    /`---'   ,' \n\
        ;    `\"\"\"\"'   (  \\:::'     /     _,-'\n\
        ;              \\  \\:'    ,';:.,-' \n\
        (              :  )\\    ( \n\
         `.             \\   \\    ; \n\
   -hrr-   `-.___       : ,\\ \\  ( \n\
              ,','._::::| \\ \\ \\  \\ \n\
             (,(,---;;;;;  \\ \\|;;;) \n\
                         `._\\_\\");

	AddPiece("                                                 ,::::.._\n\
                                               ,':::::::::.\n\
                                           _,-'`:::,::(o)::`-,.._\n\
                                        _.', ', `:::::::::;'-..__`.\n\
                                   _.-'' ' ,' ,' ,\\:::,'::-`'''\n\
                               _.-'' , ' , ,'  ' ,' `:::/\n\
                         _..-'' , ' , ' ,' , ,' ',' '/::\n\
                 _...:::'`-..'_, ' , ,'  , ' ,'' , ,'::|\n\
              _`.:::::,':::::,'::`-:..'_',_'_,'..-'::,'|\n\
      _..-:::'::,':::::::,':::,':,'::,':::,'::::::,':::;\n\
        `':,'::::::,:,':::::::::::::::::':::,'::_:::,'/\n\
        __..:'::,':::::::--''' `-:,':,':::'::-' ,':::/\n\
   _.::::::,:::.-''-`-`..'_,'. ,',  , ' , ,'  ', `','\n\
 ,::SSt:''''`                 \\:. . ,' '  ,',' '_,'\n\
                               ``::._,'_'_,',.-'\n\
                                   \\\\ \\\\\n\
                                    \\\\_\\\\\n\
                                     \\\\`-`.-'_\n\
                                  .`-.\\\\__`. ``\n\
                                     ``-.-._\n\
                                         `");
	SelectOne();
	BTOUT;
}


CAsciiArt::~CAsciiArt()
{
}


void CAsciiArt::AddPiece(const char *s)
{
	BTIN;
	CAsciiPiece *p;
	const char *q;
	int tw, th, z;

	try { p = new CAsciiPiece(); } catch(...) { p = NULL; }
	if (p == NULL) NewException((double)sizeof(CAsciiPiece),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	m_oaAsciiPieces.Add(p);
	p->m_iWidth = 0;
	tw = 0;
	p->m_iHeight = 0;
	q = s;
	while (*q != 0)
	{
		q++;
		tw++;
		if ((*q == '\n') || (*q == 0))
		{
			if (tw-1 > p->m_iWidth)
				p->m_iWidth = tw-1;
			tw = 0;
			p->m_iHeight++;
			if (*q == 0)
				break;
		}
	}

	try { p->m_pBuf = new char[p->m_iWidth*p->m_iHeight]; } catch(...) { p->m_pBuf = NULL; }
	if (p->m_pBuf == NULL) NewException((double)p->m_iWidth*p->m_iHeight*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	tw = 0;
	th = 0;
	q = s;
	while (*q != 0)
	{
		p->m_pBuf[tw+th*p->m_iWidth] = *q;
		q++;
		tw++;
		if ((*q == '\n') || (*q == 0))
		{
			for (z=tw;z<p->m_iWidth;z++)
				p->m_pBuf[z+th*p->m_iWidth] = ' ';
			tw = 0;
			th++;
			if (*q == 0)
				break;
			q++;
		}
	}
/*	printf("### Anfang ###\n");
	for (th=0;th<p->m_iHeight;th++)
	{
		printf("\"");
		for (tw=0;tw<p->m_iWidth;tw++)
			printf("%c",p->m_pBuf[tw+th*p->m_iWidth]);
		printf("\"\n");
	}
	printf("### Ende ###\n");*/
	BTOUT;
}


void CAsciiArt::SelectOne()
{
	BTIN;
	m_iSelectedPiece = rand()%m_oaAsciiPieces.GetSize();
	BTOUT;
}


int CAsciiArt::GetHeight()
{
	BTIN;
	if (m_iSelectedPiece == -1)
		SelectOne();
	BTOUT;
	return ((CAsciiPiece*)m_oaAsciiPieces[m_iSelectedPiece])->m_iHeight;
}


int CAsciiArt::GetWidth()
{
	BTIN;
	if (m_iSelectedPiece == -1)
		SelectOne();
	BTOUT;
	return ((CAsciiPiece*)m_oaAsciiPieces[m_iSelectedPiece])->m_iWidth;
}


char CAsciiArt::GetChar()
{
	BTIN;
	char c;
	if (m_iSelectedPiece == -1)
		SelectOne();
	if (m_iPosY == -1)
	{
		if (m_iPosX == -1)
			c = '/';
		else if (m_iPosX == ((CAsciiPiece*)m_oaAsciiPieces[m_iSelectedPiece])->m_iWidth)
			c = '\\';
		else if (m_iPosX > ((CAsciiPiece*)m_oaAsciiPieces[m_iSelectedPiece])->m_iWidth)
			c = 0;
		else 
			c = '=';
	} else if (m_iPosY == ((CAsciiPiece*)m_oaAsciiPieces[m_iSelectedPiece])->m_iHeight+1)
	{
		if (m_iPosX > 0)
			c = 0;
		else c = ' ';
	} else if (m_iPosY == ((CAsciiPiece*)m_oaAsciiPieces[m_iSelectedPiece])->m_iHeight)
	{
		if (m_iPosX == -1)
			c = '\\';
		else if (m_iPosX == ((CAsciiPiece*)m_oaAsciiPieces[m_iSelectedPiece])->m_iWidth)
			c = '/';
		else if (m_iPosX > ((CAsciiPiece*)m_oaAsciiPieces[m_iSelectedPiece])->m_iWidth)
			c = 0;
		else 
			c = '=';
	} else if (m_iPosX == ((CAsciiPiece*)m_oaAsciiPieces[m_iSelectedPiece])->m_iWidth)
		c = '|';
	else if (m_iPosX == -1)
		c = '|';
	else if (m_iPosX > ((CAsciiPiece*)m_oaAsciiPieces[m_iSelectedPiece])->m_iWidth)
		c = 0;
	else c = ((CAsciiPiece*)m_oaAsciiPieces[m_iSelectedPiece])->GetAt(m_iPosX,m_iPosY);
	m_iPosX++;
	BTOUT;
	return c;
}


char CAsciiPiece::GetAt(int x, int y)
{
	return m_pBuf[x+y*m_iWidth];
}


void CAsciiArt::NewLine()
{
	BTIN;
	if (m_iSelectedPiece == -1)
		SelectOne();
	m_iPosX = -1;
	m_iPosY++;
	if (m_iPosY > ((CAsciiPiece*)m_oaAsciiPieces[m_iSelectedPiece])->m_iHeight+1)
	{
		SelectOne();
		m_iPosY = -1;
	}
	BTOUT;
}
