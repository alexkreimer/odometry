%latex2e declaration
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{iitcsthesis}[2008/11/25 Technion cs thesis class]
% Written by Ido Tal, on the basis of Dan Kenigsberg's PhD thesis.
% Thanks to Gabi Kliot for a much needed fix.

\LoadClass[11pt,twoside]{report}

%\RequirePackage[hebrew,english]{babel}
\RequirePackage{etoolbox}
\RequirePackage{fontspec}
\RequirePackage{setspace}

\renewcommand{\baselinestretch}{1.1}

\newcommand{\@authorEnglish}{}
\newcommand{\@authorHebrew}{}
\newcommand{\@titleEnglish}{}
\newcommand{\@titleHebrew}{}
\newcommand{\@supervisorEnglish}{}
\newcommand{\@supervisorHebrew}{}
\newcommand{\@GregorianDateEnglish}{}
\newcommand{\@GregorianDateHebrew}{}
\newcommand{\@JewishDateEnglish}{}
\newcommand{\@JewishDateHebrew}{}
\newcommand{\@personalThanksEnglish}{}
\newcommand{\@personalThanksHebrew}{}
\newcommand{\@financialThanksEnglish}{}
\newcommand{\@financialThanksHebrew}{}

\newcommand{\authorEnglish}[1]{\renewcommand{\@authorEnglish}{#1}}
\newcommand{\authorHebrew}[1]{\renewcommand{\@authorHebrew}{#1}}
\newcommand{\titleEnglish}[1]{\renewcommand{\@titleEnglish}{#1}}
\newcommand{\titleHebrew}[1]{\renewcommand{\@titleHebrew}{#1}}
\newcommand{\supervisorEnglish}[1]{\renewcommand{\@supervisorEnglish}{#1}}
\newcommand{\supervisorHebrew}[1]{\renewcommand{\@supervisorHebrew}{#1}}
\newcommand{\GregorianDateEnglish}[1]{\renewcommand{\@GregorianDateEnglish}{#1}}
\newcommand{\GregorianDateHebrew}[1]{\renewcommand{\@GregorianDateHebrew}{#1}}
\newcommand{\JewishDateEnglish}[1]{\renewcommand{\@JewishDateEnglish}{#1}}
\newcommand{\JewishDateHebrew}[1]{\renewcommand{\@JewishDateHebrew}{#1}}
\newcommand{\personalThanksEnglish}[1]{\renewcommand{\@personalThanksEnglish}{#1}}
\newcommand{\personalThanksHebrew}[1]{\renewcommand{\@personalThanksHebrew}{#1}}
\newcommand{\financialThanksEnglish}[1]{\renewcommand{\@financialThanksEnglish}{#1}}
\newcommand{\financialThanksHebrew}[1]{\renewcommand{\@financialThanksHebrew}{#1}}

\DeclareOption{PhD}{%
\newcommand{\@degreeEnglish}{Doctor of Philosophy}
\newcommand{\@degreeHebrew}{דוקטור לפילוסופיה}
}

\DeclareOption{MSc}{%
\newcommand{\@degreeEnglish}{Master of Science in Computer Science}
\newcommand{\@degreeHebrew}{מגיסטר למדעים במדעי המחשב}
}

\DeclareOption{beforeExam}{% No personal thanks are printed before the exam is taken and passed
\renewcommand{\personalThanksEnglish}[1]{\renewcommand{\@personalThanksEnglish}{}}
\renewcommand{\personalThanksHebrew}[1]{\renewcommand{\@personalThanksHebrew}{}}
}

\newcommand{\maketitleEnglish}{%
\pagestyle{empty}
%cover page
\begin{center}
\begin{spacing}{2.5}
{\Huge \bf \@titleEnglish}
\end{spacing}
\vfill
{\LARGE \bf \@authorEnglish}
\vfill
\cleardoublepage
\end{center}

\begin{center}
\begin{spacing}{2.5}
{\Huge \bf \@titleEnglish}
\end{spacing}
\Large
\vfill
Research Thesis
\vfill
{\normalsize
Submitted in partial fulfillment of the requirements \\
for the degree of \@degreeEnglish }
\vfill
{\LARGE \bf \@authorEnglish}
\vfill
\parbox{3.5in}{
\normalsize
{\center Submitted to the Senate of \\
the Technion --- Israel Institute of Technology \\ }
\@JewishDateEnglish \hfill Haifa \hfill \@GregorianDateEnglish
}
\end{center}

\cleardoublepage
%acknowledgment page
\noindent
\@supervisorEnglish
\vfill
\noindent
\@personalThanksEnglish
\vfill
\noindent
\@financialThanksEnglish
\cleardoublepage

\pagestyle{plain}
\pagenumbering{roman}

\tableofcontents
\listoffigures
%\listoftables
\cleardoublepage
}

\newcommand{\abstractEnglish}{%
\pagenumbering{arabic}

\chapter*{Abstract}
\addcontentsline{toc}{chapter}{Abstract}
\addcontentsline{cot}{chapter}{תקציר}
\setcounter{page}{1}
}

\newcommand{\abbreviationsAndNotationsEnglish}{%
\chapter*{Abbreviations and Notations}
\addcontentsline{toc}{chapter}{Abbreviations and Notations}
\addcontentsline{cot}{chapter}{רשימת קיצורים וסמלים}
}

\newcommand{\maketitleHebrew}{
\cleardoublepage

\begin{otherlanguage}{hebrew}
\pagestyle{empty}

%cover page
\begin{center}
\begin{spacing}{2.5}
{\Huge \bf \@titleHebrew}
\end{spacing}
\vfill
{\LARGE \bf \@authorHebrew}
\vfill
\cleardoublepage
\end{center}

\begin{center}
\Large
\begin{spacing}{2.5}
{\Huge \bf \@titleHebrew}
\end{spacing}
\vfill
חיבור על מחקר
\vfill
{\normalsize
לשם מילוי חלקי של הדרישות לקבלת התואר \\
\@degreeHebrew}
\vfill
{\LARGE \bf \@authorHebrew}
\vfill
\parbox{3in}{
\normalsize
{\center הוגש לסנט הטכניון \L{--} מכון טכנולוגי לישראל \\ }
\@JewishDateHebrew
\hfill חיפה \hfill
\@GregorianDateHebrew
}
\end{center}

\cleardoublepage

%acknowledgment page
\noindent
\@supervisorHebrew
\vfill
\noindent
\@personalThanksHebrew
\vfill
\noindent
\@financialThanksHebrew
\cleardoublepage

\pagestyle{plain}
\pagenumbering{alph}

\cleardoublepage
}


\newcommand{\abstractHebrew}{%
\chapter*{תקציר}
\addcontentsline{cot}{chapter}{תקציר בעברית}
\addcontentsline{toc}{chapter}{Abstract in Hebrew}
\setlength{\parindent}{0pt}
\setlength{\parskip}{10pt}
}

\newcommand{\finishedHebrew}{%
\clearpage

\end{otherlanguage}
}

\ProcessOptions\relax
