%!TeX root = rapportMaster.tex

\NeedsTeXFormat{LaTeX2e}
\ProvidesPackage{isddmaster}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
\LoadClass[11pt,twoside]{article}


\def\mymaster{BI-ISDD}


%---- Spacing around titles
\renewcommand\section{\@startsection {section}{1}{\z@}%
	{-5.0ex \@plus -1ex \@minus -.2ex}%  before skip
	{1.3ex \@plus.2ex}% after skip
	{\normalfont\Large\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
	{-3.0ex\@plus -1ex \@minus -.2ex}% 
	{0.5ex \@plus .2ex}%
	{\normalfont\large\bfseries}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
	{-2.0ex\@plus -1ex \@minus -.2ex}% 
	{0.3ex \@plus .2ex}%
	{\normalfont\normalsize\bfseries}}

%---- Table of Contents and sectioning depth
\setcounter{secnumdepth}{3} 
\setcounter{tocdepth}{2}


%---- Core packages (each loaded only once)
\RequirePackage{graphicx}
\RequirePackage[x11names]{xcolor}

\RequirePackage[pdfstartview=FitH,%
bookmarks=true,%
bookmarksopen=true,%
breaklinks=true,%
colorlinks=true,%
linkcolor=black,anchorcolor=blue,%
citecolor=black,filecolor=blue,%
menucolor=blue,%
urlcolor=blue]{hyperref}

\RequirePackage{booktabs}
\RequirePackage[scaled=.90]{helvet}

% Font stack: use lmodern as safe fallback (always installed)
\RequirePackage{lmodern}
\RequirePackage[T1]{fontenc}

\RequirePackage{amsmath}
\RequirePackage{lineno}
\RequirePackage[a4paper,margin=2.5cm]{geometry}

\parskip 5em plus1pt minus1pt
\parindent 5em

%---- Bibliography: bibtex/natbib only (compatible with \bibliographystyle + \bibliography)
\RequirePackage[numbers,sort&compress]{natbib}

%---- Other packages
\RequirePackage{caption}
\RequirePackage{subfig}
\RequirePackage{multirow}
\RequirePackage{multicol}
\RequirePackage{wrapfig}
\RequirePackage{framed}
\RequirePackage{setspace}
\RequirePackage{emptypage}

%---- Listings
\RequirePackage{listings}
\RequirePackage{etoolbox}

\definecolor{backcolour}{rgb}{0.95,0.95,0.92}
\lstdefinestyle{mystyle}{
	backgroundcolor=\color{backcolour},
	basicstyle=\ttfamily\footnotesize,
	breakatwhitespace=false,
	breaklines=true,
	captionpos=b,
	keepspaces=true,
	numbers=left,
	numbersep=5pt,
	numberstyle=\footnotesize\color{black},
	showspaces=false,
	showstringspaces=false,
	showtabs=false,
	tabsize=2,
	commentstyle=\color{black},
	keywordstyle=\color{black},
	stringstyle=\color{black},
	identifierstyle=\color{black},
}
\lstset{style=mystyle, language=bash}

\setlength{\columnsep}{10pt}

%---- Title page defaults
\def\author{Anna Perova}
\def\title{La dynamique de la protease VIH-2 (PR2): \\[5mm] Rapport de Stage M1 }
\def\shorttitle{Rapport de stage}
\def\unit{BFA CNRS UMR 8251 - \\ INSERM ERL 1133}
\def\team{Profilage Pharmacologique In Silico (IsPP)}
\def\spe{Modélisation des macromolécules}
\def\supervisor{L. Regad \& M. Baillif}
\def\mydate{\today}
\def\keywords{Sciences de la Terre, master, rapport, conseils}

\newcounter{firstyear}
\newcounter{lastyear}[firstyear]
\newcounter{thismonth}

\newcommand{\checkyears}{
	\setcounter{firstyear}{\the\year}
	\setcounter{thismonth}{\the\month}
	\ifnum\the\value{thismonth}<9
	\addtocounter{firstyear}{-1}
	\fi
	\setcounter{lastyear}{\value{firstyear}}
	\stepcounter{lastyear}
}

\newcommand{\yearlabel}{Année Universitaire }
\newcommand{\unitlabel}{Unité de Recherche : }
\newcommand{\teamslabel}{\'{E}quipe(s) : }
\newcommand{\suplabel}{Encadrant(s) : }
\newcommand{\abstractlabel}{Résumé}
\newcommand{\keyslabel}{Mots clés : }

%---- French/English
\RequirePackage[english,french]{babel}
\RequirePackage{iflang}
\def\frenchcontentsname{Sommaire}

\newcommand{\Entete}{
	
	\IfLanguagePatterns{english}{
		\renewcommand{\yearlabel}{Academic Year }
		\renewcommand{\unitlabel}{Research Unit: }
		\renewcommand{\teamslabel}{Team(s): }
		\renewcommand{\suplabel}{Supervisor(s): }
		\renewcommand{\abstractlabel}{Abstract}
		\renewcommand{\keyslabel}{Keywords: }
	}{}
	
	\vspace*{-2cm}
	
	\begin{minipage}{0.6\linewidth}
		\yearlabel \thefirstyear\, -- \thelastyear\\
		\unitlabel  \unit\\
		\teamslabel \team\\
		\suplabel  \supervisor\\
	\end{minipage}
	\hfill
	\medskip
	\begin{minipage}{0.30\linewidth}
		\includegraphics[width=\linewidth]{Logo/logoUPC.png}\\[2mm]
		\includegraphics[width=\linewidth]{Logo/mymaster.png}\\[2mm]
		\hspace*{\fill}{\sc \spe}\hspace*{1.5mm}\\[0.2cm]
	\end{minipage}
	\smallskip
	{\parskip 3pt
		\begin{center}
			\hrule
			\vspace*{1mm}
			{\huge\bfseries \title}\\[5mm]
			\hrule
			\vspace*{4mm}
			{\Large\bfseries \author}\\[3mm]
			\mydate
		\end{center}
}}

%---- Headers and footers
\RequirePackage{fancyhdr}
\addtolength{\headheight}{1.6pt}
\lhead[\fancyplain{}{\bfseries\thepage}]{\fancyplain{}{\author}}
\rhead[\fancyplain{}{\nouppercase{\bfseries \shorttitle}}]%
{\fancyplain{}{\bfseries\thepage}}
\fancyfoot{}
\fancyfoot[LE,RO]{\itshape \yearlabel \thefirstyear\, -- \thelastyear}
\pagestyle{fancyplain}

%---- Abstract and Keywords
\renewenvironment{abstract}
{
	\thispagestyle{empty}
	\begin{center}
		{\bfseries \abstractlabel}
	\end{center}
	\small
}
{
	\par
	{\bfseries \keyslabel} \keywords
	\par}



%---- Custom captions (Figure X -- caption style)

\captionsetup{font=small}

\long\def\@caption#1[#2]#3{\par\addcontentsline{\csname
		ext@#1\endcsname}{#1}{\protect\numberline{\csname
			the#1\endcsname}{\ignorespaces #2}}\begingroup
	\@parboxrestore
	\normalsize
\ifx\@tempa\relax\let\@tempa\@makecaption\fi
	\@tempa{\csname fnum@#1\endcsname}{\ignorespaces #3}\\
	\endgroup}

\long\def\@makefigurecaption#1#2{%
	\vskip 1pt 
	\setbox\@tempboxa\hbox{#1 -- {\small #2}}%
	\ifdim \wd\@tempboxa >\hsize #1 -- {\small #2}\\ \else
	\hbox to\hsize{\hfil\box\@tempboxa\hfil}%
	\fi
}	\expandafter\let\expandafter\@tempa\csname @make#1caption\endcsname
	

\long\def\@maketablecaption#1#2{%
	\setbox\@tempboxa\hbox{#1 -- {\small #2}}%
	\ifdim \wd\@tempboxa >\hsize #1 -- {\small #2}\\ \else
	\hbox to\hsize{\hfil\box\@tempboxa\hfil}%
	\fi
	\vskip 1pt 
}