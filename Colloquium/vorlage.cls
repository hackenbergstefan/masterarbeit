%% LaTeX Beamer presentation template (requires beamer package)
%% see http://bitbucket.org/rivanvx/beamer/wiki/Home
%% idea contributed by H. Turgut Uyar
%% template based on a template by Till Tantau
%% this template is still evolving - it might differ in future releases!

\LoadClass[9pt]{beamer}
% 
\mode<presentation>
%\mode<handout>

\RequirePackage{import} %% organize imports

\setbeamertemplate{navigation symbols}{}%remove navigation symbols

\usetheme{CambridgeUS}
\hypersetup{pdfpagemode=FullScreen}

\RequirePackage[ngerman]{babel}
\RequirePackage[utf8]{inputenc}
\RequirePackage[T1]{fontenc}
\RequirePackage{amsmath,amsfonts,amssymb,wasysym,array}
\RequirePackage[sans]{dsfont}
\RequirePackage{MnSymbol} % downmapsto
\RequirePackage{rotating,xspace}
\RequirePackage{xifthen,printlen,etextools}
%\usepackage[backend=biber,style=verbose]{biblatex}
%\bibliography{biblio.bib}
%\nocite{*}
%\RequirePackage{pgfpages}
%\pgfpagesuselayout{4 on 1}[a4paper,landscape,border shrink=5mm]

\usepackage{listings}

\RequirePackage{tikz}
\usetikzlibrary{arrows,calc, positioning, fit, shapes, backgrounds, shadows,
  matrix,decorations.pathreplacing}

% \usepackage{pgfpages}
% \pgfpagesuselayout{4 on 1}[a4paper,landscape,border shrink=5mm]

% font definitions, try \usepackage{ae} instead of the following
% three lines if you don't like this look
% \usepackage{mathptmx}
% \usepackage[scaled=.90]{helvet}
% \usepackage{courier}


% \setbeameroption{show notes on second screen=left}

\renewcommand*{\thefootnote}{\arabic{footnote}}

%\AtBeginDocument{
%\setlength{\abovedisplayshortskip}{0pt}
%\setlength{\belowdisplayshortskip}{0pt}
%\setlength{\abovedisplayskip}{5pt}
%\setlength{\belowdisplayskip}{0pt}
%}



%% Define colors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\tikzstyle{emptynode}=[inner sep=0pt, outer sep=0pt]
\tikzstyle{colored}=[fill=blue!20, rounded corners=2pt, inner sep=3pt, outer sep=-3pt]
\tikzstyle{overlaynode}=[text=blue!60, fill=gray, draw=blue!20,
	rounded corners=5pt, fill opacity=0.2, text opacity=1]
\tikzstyle{overlayarrow}=[draw=blue!70, opacity=0.3, line width=1pt]

\tikzstyle{overlayproofnode}=[text=blue!60, draw=blue!20, line width=0.5pt,
	font=\scriptsize]
\tikzstyle{overlayproofarrow}=[draw=blue!60, opacity=0.5, line width=1pt,
	shorten >= 3pt]

\newcommand{\emptynode}[3][]{%
	\tikz[remember picture, baseline, outer sep=0pt, inner sep=0pt]{
		\node[anchor=base, emptynode, #1] (#2) {#3};}}
		

%% Math abkürzungen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
\newenvironment{notation}
	{\begin{exampleblock}{Notation}}	
	{\end{exampleblock}}

\newenvironment{beispiel}[1][]
	{\begin{exampleblock}{Beispiel #1}}	
	{\end{exampleblock}}

\newenvironment{satz}[1][]%
  {\begin{block}{Satz #1}}
  {\end{block}}
 
\renewenvironment{lemma}[1][]%
  {\begin{block}{Lemma #1}}
  {\end{block}}


\newcommand{\ds}[1]{{\displaystyle #1}}
\let\fr\mathfrak

\newcommand{\A}{\ensuremath \mathds{A}}
\newcommand{\R}{\ensuremath \mathcal{R}}
\newcommand{\N}{\ensuremath \mathds{N}}
\newcommand{\Q}{\ensuremath \mathds{Q}}
\newcommand{\Z}{\ensuremath \mathds{Z}}
\newcommand{\C}{\ensuremath \mathcal{C}}
\newcommand{\CN}{\ensuremath \mathcal{CN}}
\newcommand{\PN}{\ensuremath \mathcal{PN}}
\newcommand{\PCN}{\ensuremath \mathcal{PCN}}
\newcommand{\G}{\ensuremath \mathcal{G}}
\newcommand{\F}{\ensuremath \mathds{F}}
\newcommand{\K}{\ensuremath \mathds{K}}
\renewcommand{\P}{\ensuremath \mathds{P}}
\newcommand{\Kb}{\ensuremath \overline K}
\renewcommand{\O}{\ensuremath \mathcal{O}}
\renewcommand{\L}{\ensuremath \mathcal{L}}
\renewcommand{\l}{\ensuremath \ell}
\renewcommand{\S}{\ensuremath\mathcal{S}}
\newcommand{\m}{\ensuremath \mathfrak{m}}
\newcommand{\speq}[1]{\ #1\ }
\newcommand{\const}{\ensuremath \mathrm{const}}
\newcommand{\divp}[1]{\ensuremath [#1]} %Divisorpunkt
\newcommand{\probn}[1]{{\sffamily #1}} %Problem Name
\newcommand{\obda}{{\small oBdA\space}}
\newcommand{\inv}{^{-1}}
\newcommand{\kl}[1]{^{(#1)}}
\newcommand{\sage}{\texttt{Sage}\xspace}
\newcommand{\cython}{\texttt{Cython}\xspace}
\newcommand{\python}{\texttt{Python}\xspace}
\newcommand{\Clang}{\texttt{C}\xspace}
\let\cal\mathcal
\newcommand{\lbr}{[}
\newcommand{\rbr}{]}

\let\div\undefined
\DeclareMathOperator{\charak}{char}
\DeclareMathOperator{\div}{div}
\DeclareMathOperator{\ord}{ord}
\DeclareMathOperator{\summ}{sum}
\DeclareMathOperator{\comp}{\circ}
\DeclareMathOperator{\ggT}{ggT}
\DeclareMathOperator{\kgV}{kgV}
\DeclareMathOperator{\supp}{supp}
\DeclareMathOperator{\Div0}{Div^0}
\DeclareMathOperator{\Divv}{Div}
\DeclareMathOperator{\Pic0}{Pic^0}
\DeclareMathOperator{\Gal}{Gal}
\DeclareMathOperator{\End}{End}
\DeclareMathOperator{\Ord}{Ord}
\DeclareMathOperator{\im}{im}
\DeclareMathOperator{\id}{id}
\DeclareMathOperator{\Tr}{Tr}
\DeclareMathOperator{\Nm}{Nm}
\DeclareMathOperator{\spann}{span}
\DeclareMathOperator{\ideal}{\triangleleft}
\DeclareMathOperator{\Ann}{Ann}
\DeclareMathOperator{\cl}{cl}


\newcommand{\funcdef}[1]{%
    \begin{array}[t]{>{\displaystyle}r>{\displaystyle}c>{\displaystyle}l}
    #1\end{array}}
\let\xto\xrightarrow

\newcommand{\partto}{\ \tikz[baseline=-0.6ex]{
    \draw[->]
        (0,0) -- +(0.4em,0) +(0.7em,0) -- +(1.3em,0);}\ }


\let\mpar\marginpar
\renewcommand{\marginpar}[1]{
    \mpar{\normalfont\sffamily\scriptsize #1}}


\newcommand{\secframe}{%
  \begin{frame}
    \begin{tikzpicture}[remember picture,overlay]
      \node[font=\huge\bfseries, text width=0.7\textwidth, text=col1,
        anchor=west]
        (sectitle)
        at ($(current page.center)-(0.3\textwidth,0)$) {\insertsection};
      \node[anchor=east, left=1cm of sectitle.west,
        font=\huge\bfseries, text=gray] 
        {\thesection};
      \draw[line width=2pt,gray!50] 
        ($(sectitle.north west)-(0.5cm,0)$) -- 
        ($(sectitle.south west)-(0.5cm,0)$);
    \end{tikzpicture}
  \end{frame}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LstListings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\lstset{
  mathescape = true,
  basicstyle = \small\normalfont\sffamily,
  frame = tb,
  framexleftmargin = 15pt,
% numbers = left,
  numberstyle = \tiny,
% numbersep = 5pt,
  breaklines = true,
  xleftmargin = 0.1\linewidth,
  xrightmargin = 0.1\linewidth,
  escapeinside = {(*}{*)},
  tabsize=3,
  morekeywords={if, and, or, is, then, else, endif, while, endwhile, for, from,
  to, do, endfor, Input, Output, Algorithmus, return},
  morecomment=[l]{//},
  columns=flexible
}

\lstdefinestyle{mycstyle}{%
  commentstyle=\color{gray}\small\itshape,
  keywordstyle={\color{col1}},
  identifierstyle={\color{col1!50!black}}
}

\lstnewenvironment{cexample}[1][]{%
  \lstset{ basicstyle = \small\normalfont\ttfamily,
  frame = none,
  %numbers = left,
  %numberstyle = \tiny,
  %numbersep = 5pt,
  breaklines = true,
  xleftmargin = 3ex,
  xrightmargin = 0pt,
  escapeinside = ||,
  tabsize=3,
  language=C,
  style=mycstyle,
  %mathescape=true,
  #1}}{}


\lstnewenvironment{sageexample}[1][]{%
  \lstset{ basicstyle = \small\normalfont\ttfamily,
  frame = none,
  %numbers = left,
  %numberstyle = \tiny,
  %numbersep = 5pt,
  breaklines = true,
  xleftmargin = 3ex,
  xrightmargin = 0pt,
  escapeinside=||,
  tabsize=3,
  language=Python,
  commentstyle=\color{gray}\small\itshape,
  %mathescape=true,
  #1}}{}



\newcommand{\ttgray}{\color{gray}\ttfamily}

\newcommand{\lstnospace}{\vspace{-\baselineskip}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Template %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\colorlet{col1}{red!60!black}
\colorlet{col2}{blue!30}

\setbeamertemplate{enumerate item}{%
  \tikz[baseline,gray]{
    \node[font=\scriptsize, anchor=base] (a)
      {\insertenumlabel}; 
    \begin{scope}[overlay]
      \draw[line width=1pt,gray!20,xshift=2pt] (a.north east) -- (a.south east);
    \end{scope}
  }}

\setbeamertemplate{itemize item}{%
  \tikz[baseline=-2pt]{
    \fill[gray] (0,0) circle (1pt);
  }}
\setbeamertemplate{itemize subitem}{%
  \tikz[baseline=-2pt]{
    \fill[gray] (0,0) circle (0.8pt);
  }}

\setbeamercolor{section in toc}{fg=black,bg=white}
\setbeamercolor{alerted text}{fg=darkred!80!gray}
\setbeamercolor*{palette primary}{fg=col1,bg=gray!30!white}
\setbeamercolor*{palette secondary}{fg=col1,bg=gray!15!white}
\setbeamercolor*{palette tertiary}{bg=col1,fg=gray!10!white}
\setbeamercolor*{palette quaternary}{fg=darkred,bg=gray!5!white}

\setbeamercolor*{sidebar}{fg=darkred,bg=gray!15!white}

\setbeamercolor*{palette sidebar primary}{fg=col1}
\setbeamercolor*{palette sidebar secondary}{fg=white}
\setbeamercolor*{palette sidebar tertiary}{fg=col1}
\setbeamercolor*{palette sidebar quaternary}{fg=gray!10!white}

%\setbeamercolor*{titlelike}{parent=palette primary}
\setbeamercolor{titlelike}{parent=palette primary,fg=darkred}
\setbeamercolor{frametitle}{bg=gray!10!white}
\setbeamercolor{frametitle right}{bg=gray!60!white}

\setbeamercolor*{separation line}{}
\setbeamercolor*{fine separation line}{}

