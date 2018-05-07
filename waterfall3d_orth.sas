/******************************************************************************
*  Program name     : waterfall3d.sas
*  Project          : 
*  Written by       : Liz Thomas
*  Date of creation : 09Aug2017
*  Description      : Create 3d Waterfall Plot
*  Output data      :
*  Macro call       :
*  Revision History :
*  Date:        Author          Description of the change
*******************************************************************************/


/* Project 3d coordinates into 2d space for plotting
  Based off of code here: http://blogs.sas.com/content/graphicallyspeaking/2015/03/10/a-3d-scatter-plot-macro/
*/

/* Functions from sourced link: for matrix inversion, matrix multiplication and creation of identity matrix */
options cmplib=sasuser.funcs;

proc fcmp outlib=sasuser.funcs.mat;
  subroutine MatInv(Mat[*,*], InvMat[*,*]);
  outargs InvMat;
  call inv(Mat, InvMat);
  endsub;

  subroutine MatMult(A[*,*], B[*,*], C[*,*]);
  outargs C;
  call mult(A, B, C);
  endsub;

  subroutine MatIdent(A[*,*]);
  outargs A;
  call identity(A);
  endsub;
run;
quit;

/* Macro to normalize 3d coordinate dataset to unit cube plot -> maps range to [-1,1]. */
%macro normalize(dataset=, x=x, y=y, z=z, recalclims=1);
%if &recalclims=1 %then %do;
proc sql;
	create table lims as
	select min(x) as xmin, max(x) as xmax, min(y) as ymin, max(y) as ymax, min(z) as zmin, max(z) as zmax
	from &dataset.
	;
quit;
data _null_;
	set lims;
	if (missing(&xmin.) or xmin<&xmin.) then call symput("xmin",floor(xmin));
	if (missing(&xmax.) or xmax>&xmax.) then call symput("xmax",ceil(xmax));
	if (missing(&ymin.) or ymin<&ymin.) then call symput("ymin",floor(ymin));
	if (missing(&ymax.) or ymax>&ymax.) then call symput("ymax",ceil(ymax));
	if (missing(&zmin.) or zmin<&zmin.) then call symput("zmin",floor(zmin));
	if (missing(&zmax.) or zmax>&zmax.) then call symput("zmax",ceil(zmax));
run;
%end;
data &dataset._n(drop=xrange yrange zrange);
	set &dataset;
	xrange=&xmax-&xmin;
	yrange=&ymax-&ymin;
	zrange=&zmax-&zmin;
	x=0.99*(2*(&x-&xmin)/xrange -1)+0.005;
	y=0.99*(2*(&y-&ymin)/yrange-1)+0.005;
	z=0.99*(2*(&z-&zmin)/zrange-1)+0.005;
run;
%mend;


/* Projection macro -> translates &x, &y, &z into coordinates x2 and y2 */
%macro project(ds=, x=x, y=y, z=z, rotx=100, roty=200, rotz=180);
data &ds._p;
  array u[4,4] _temporary_;  /*--Intermediate Matrix--*/
  array v[4,4] _temporary_;  /*--Intermediate Matrix--*/
  array uu[4,4] _temporary_; /*--Intermediate Matrix--*/
  array w[4,4] _temporary_;  /*--Final View Matrix--*/
  array m[4,4] _temporary_;  /*--Projection Matrix--*/
  array rx[4,4] _temporary_; /*--X rotation Matrix--*/
  array ry[4,4] _temporary_; /*--Y rotation Matrix--*/
  array rz[4,4] _temporary_; /*--Z rotation Matrix--*/
  array tr[4,4] _temporary_; /*--Translation Matrix--*/
  array d[4,1] _temporary_;  /*--World Data Array --*/
  array p[4,1] _temporary_;  /*--Projected Data Array --*/
  retain r t f n;
  r=1; t=1; f=1; n=-1;
  pi=constant("PI");
  fac=pi/180;
  A=&rotx.*fac; B=&roty.*fac; C=&rotz.*fac;

  /*--Set up orthographic projection matrix--*/
  m[1,1]=1/r;   m[1,2]=0.0;  m[1,3]=0.0;      m[1,4]=0.0;
  m[2,1]=0.0;   m[2,2]=1/t;  m[2,3]=0.0;      m[2,4]=0.0;
  m[3,1]=0.0;   m[3,2]=0.0;  m[3,3]=-2/(f-n); m[3,4]=-(f+n)/(f-n);
  m[4,1]=0.0;   m[4,2]=0.0;  m[4,3]=0.0;      m[4,4]=1.0;

  /*--Set up translation matrix--*/
  tr[1,1]=1.0;   tr[1,2]=0.0;  tr[1,3]=0.0;      tr[1,4]=0.0;
  tr[2,1]=0.0;   tr[2,2]=1.0;  tr[2,3]=0.0;      tr[2,4]=0.0;
  tr[3,1]=0.0;   tr[3,2]=0.0;  tr[3,3]=1.0;		 tr[3,4]=0.0;
  tr[4,1]=0.0;   tr[4,2]=0.0;  tr[4,3]=0.0;      tr[4,4]=1.0; 

  /*--Set up X rotation matrix--*/
  rx[1,1]=1;     rx[1,2]=0.0;     rx[1,3]=0.0;      rx[1,4]=0.0;
  rx[2,1]=0.0;   rx[2,2]=cos(A);  rx[2,3]=-sin(A);  rx[2,4]=0.0;
  rx[3,1]=0.0;   rx[3,2]=sin(A);  rx[3,3]=cos(A);   rx[3,4]=0.0;
  rx[4,1]=0.0;   rx[4,2]=0.0;     rx[4,3]=0.0;      rx[4,4]=1.0;

  /*--Set up Y rotation matrix--*/
  ry[1,1]=cos(B);  ry[1,2]=0.0;  ry[1,3]=sin(B);  ry[1,4]=0.0;
  ry[2,1]=0.0;     ry[2,2]=1.0;  ry[2,3]=0.0;     ry[2,4]=0.0;
  ry[3,1]=-sin(B); ry[3,2]=0.0;  ry[3,3]=cos(B);  ry[3,4]=0.0;
  ry[4,1]=0.0;     ry[4,2]=0.0;  ry[4,3]=0.0;     ry[4,4]=1.0;

  /*--Set up Z rotation matrix--*/
  rz[1,1]=cos(C);  rz[1,2]=-sin(C); rz[1,3]=0.0;  rz[1,4]=0.0;
  rz[2,1]=sin(C);  rz[2,2]=cos(C);  rz[2,3]=0.0;  rz[2,4]=0.0;
  rz[3,1]=0.0;     rz[3,2]=0.0;     rz[3,3]=1.0;  rz[3,4]=0.0;
  rz[4,1]=0.0;     rz[4,2]=0.0;     rz[4,3]=0.0;  rz[4,4]=1.0;
  
  /*--Build transform matrix--*/
  call MatMult(ry,rx,u);    *Rotate in X direction first, then Y direction;
  call MatMult(rz, u, uu); *Rotate in Z direction;
  call MatMult(tr, uu, v);  *Translate away from the Z-origin so viewpoint is not within the field;
  call MatMult(m, v, w);   *Apply the projection matrix;

  set &ds.;

  /*--Transform data--*/
  d[1,1]=&x; d[2,1]=&y; d[3,1]=&z; d[4,1]=1;
  call MatMult(w, d, p);
  x2=p[1,1]/p[4,1]; y2=p[2,1]/p[4,1]; /*z2=p[3,1]; w2=p[4,1]; */
  rv1=p[1,1]; rv2=p[2,1]; rv3=p[3,1]; rv4=p[4,1];
run;
%mend;

/* Create 'realistic' fake data */
data mydat;
	call streaminit(59483);
	do subjid=1 to 20;
		uni1=rand('UNIFORM');
		uni2=rand('UNIFORM');
		time=5*(3*(1-uni1)+uni2);
		pctchg=max(-100,150*(uni1-0.75));
		if(uni1>0.75) then pctchg=400*(uni1-0.75);
		else pctchg=max(-100,150*(uni1-0.75));
		if abs(pctchg)<0.5 then pctchg=pctchg+1.9;
		cont=rand('BERNOULLI',0.25);
		doselvl=ceil(uni2*2);
		output;
	end;
run;

/* Cohort format */
proc format;
	value cht
	1 = "Cohort 1"
	2 = "Cohort 2"
	;
run;

proc sort data=mydat;
	by descending pctchg;
run;

proc sort data=mydat;
	by descending pctchg;
run;

data mydat;
	set mydat;
	subjid=_n_;
run;

%let delta=0.35;
%let wide=1;
%let x=subjid;
%let y=time;
%let z=pctchg;
%let xmin=.;
%let xmax=.;
%let ymin=0;
%let ymax=.;
%let zmin=-100;
%let zmax=100;
%let yarrowmult=0.05;
%let zarrowmult=0.1;

proc sql;
	create table lims as
	select min(&x.) as xmin, max(&x.) as xmax, min(&y.) as ymin, max(&y.) as ymax, min(&z.) as zmin, max(&z.) as zmax
	from mydat
	;
quit;
data _null_;
	set lims;
	if (missing(&xmin.)) then call symput("xmin",floor(xmin));
	if (missing(&xmax.)) then call symput("xmax",ceil(xmax));
	if (missing(&ymin.)) then call symput("ymin",floor(ymin));
	if (missing(&ymax.)) then call symput("ymax",ceil(ymax));
	if (missing(&zmin.)) then call symput("zmin",floor(zmin));
	if (missing(&zmax.)) then call symput("zmax",ceil(zmax));
run;

/* Polygon data in 3-D coordinates */
data polydata(drop=delta2);
	set mydat;
	length grpid delta2 8;
	delta2=&wide.*(0.5) + (1-&wide.)*&delta.;
	grpid=1;
	x=&x.-&delta; y=0; z=0; output;
	y=&y.; output;
	if cont=1 then do;
		x=&x.-delta2; output;
		x=&x.; y=&y.+&yarrowmult.*&ymax.; output;
		x=&x.+delta2; y=&y.; output;
	end;
	x=&x.+&delta; output;
	y=0; output;
	x=&x.-&delta; output;
	grpid=2;
	x=&x.-&delta; y=0; z=0; output;
	z=min(&z.,&zmax.); output;
	if &z.>&zmax. then do;
		x=&x.-delta2; output;
		x=&x.; z=(1+&zarrowmult.)*&zmax.; output;
		x=&x.+delta2; z=&zmax.; output;
	end;
	x=&x.+&delta; output;
	z=0; output;
	x=&x.-&delta; output;
run;

%normalize(dataset=polydata);

/* Axis data in 3-D coordinates */
data axes;
	length idlab x y z 8 type $5;
	type="start"; idlab=1; x=&xmin.; y=0; z=0; output;
	idlab=2; x=0; y=&ymin.; z=0; output;
	idlab=3; x=0; y=0; z=&zmin.; output;
	idlab=4; x=&xmax.; y=0; z=0; output;
	type="end"; idlab=1; x=&xmax.; y=0; z=0; output;
	idlab=2; x=0; y=&ymax.; z=0; output;
	idlab=3; x=0; y=0; z=&zmax.; output;
	idlab=4; x=&xmax.; y=&ymax.; z=0; output;
run;


/* Tickmark data */
data zticks;
	input idlab z;
	zorig=z;
	datalines;
1 -100
2 -50
3 50
4 100
;
run;

data yticks(drop=i);
	idlab=1;
	do i=1 to ceil(&ymax.)/5;
		y=5*i;
		yorig=y; 
		output; 
		idlab=idlab+1;
	end;
run;

data xtickleft;
	 type="start"; x=(&xmin.-(0.025*(&xmax.-&xmin.))); output;
	 type="end"; x=&xmin.; output;
run;

data xtickright;
	type="start"; x=(&xmax.+(0.025*(&xmax.-&xmin.))); output;
	type="end"; x=&xmax.; output;
run;

data xmid;
	type="start"; x=&xmin.; mid=1; output;
	type="end"; x=&xmax.; mid=1; output;
run;

proc sql;
	create table tick1 as
	select *, 0 as y, 1 as id2 
	from xtickleft, zticks
	;
	create table tick2 as
	select *, 0 as y, 2 as id2
	from xmid, zticks
	;
	create table tick3 as
	select *, 0 as z, 3 as id2
	from xtickright, yticks
	;
	create table tick4 as
	select *, 0 as z, 4 as id2
	from xmid, yticks
	;
quit;

data tickmarks;
	set tick1 tick2 tick3 tick4;
run;

%normalize(dataset=axes, recalclims=0);
%normalize(dataset=tickmarks, recalclims=0);

%project(ds=polydata_n);
%project(ds=axes_n);
%project(ds=tickmarks_n);

/* Annotation dataset for text labels */
data anno(drop=x2 y2);
	length label $50 function x1space y1space $20;
	set tickmarks_n_p;
	where id2 in (1,3);
	function="text";
	x1space="datavalue";
	y1space="datavalue";
	if id2=1 then do;
		x1=x2-0.10;
		y1=y2;
		label=put(zorig,best5.);
		if type='end' then delete;
	end;
	if id2=3 then do;
		x1=x2+0.05;
		y1=y2;
		label=put(yorig,best5.);
		if type='start' then delete;
	end;
run;

data anno2;
	set anno end=eof;
	output;
	if eof then do;
		x1space="graphpercent";
		y1space="graphpercent";
		function="text";
		x1=15;
		y1=50;
		label="Percent Change from Baseline";
		rotate=93;
		width=50;
		output;
		x1=80;
		y1=60;
		label="Time on study";
		rotate=40;
		output;
	end;
run;


data poly_y(rename=(x2=xy y2=yy));
	set polydata_n_p;
	where grpid=1;
run;

data poly_z(rename=(x2=xz y2=yz));
	set polydata_n_p;
	where grpid=2;
run;

data axisstart(keep=idlab xs ys);
	set axes_n_p(rename=(x2=xs y2=ys));
	where type="start";
run;

data axisend(keep=idlab xe ye);
	set axes_n_p(rename=(x2=xe y2=ye));
	where type="end";
run;

/* Solid ticks */
data tickstart1(keep=id2 idlab xs ys);
	set tickmarks_n_p(rename=(x2=xs y2=ys));
	where type="start" and mid<1;
run;

data tickend1(keep=id2 idlab xe ye);
	set tickmarks_n_p(rename=(x2=xe y2=ye));
	where type="end" and mid<1;
run;

/* Dotted reference lines */
data tickstart2(keep=id2 idlab xs2 ys2);
	set tickmarks_n_p(rename=(x2=xs2 y2=ys2));
	where type="start" and mid=1;
run;

data tickend2(keep=id2 idlab xe2 ye2);
	set tickmarks_n_p(rename=(x2=xe2 y2=ye2));
	where type="end" and mid=1;
run;

data tick_c1;
	merge tickstart1 tickend1;
	by id2 idlab;
run;

data tick_c2;
	merge tickstart2 tickend2;
	by id2 idlab;
run;

data axes_c;
	merge axisstart axisend;
	by idlab;
run;

data combined;
	set poly_y poly_z axes_c tick_c1 tick_c2;
	retain dl2;
	if doselvl>. then dl2=doselvl;
	if dl2=. then dl2=1;
run;

data dlord;
	do dl2=1 to 2;
	output;
	end;
run;

data combined;
	set dlord combined;
	label dl2 ="Dose level";
run;
	

*ods pdf file="wtfl3d_orthag.pdf";
proc sgplot data=combined nowall noborder aspect=1 noautolegend nocycleattrs sganno=anno2;
	vector x=xe y=ye / xorigin=xs yorigin=ys noarrowheads lineattrs=(color=black) attrid=AxisTick;
	vector x=xe2 y=ye2 / xorigin=xs2 yorigin=ys2 noarrowheads lineattrs=(pattern=dot color=black) attrid=RefLine;
	polygon id=&x. x=xy y=yy / fill fillattrs=(color=gold);
	polygon id=&x. x=xy y=yy / lineattrs=(color=black pattern=solid);
	polygon id=&x. x=xz y=yz / fill group=dl2 name="dose";
	polygon id=&x. x=xz y=yz / lineattrs=(color=black pattern=solid);
	keylegend "dose"/ title="Cohort";
	xaxis display=none;
	yaxis display=none;
	styleattrs datacolors=(CX58278A CX979AA7);
run;
*ods pdf close;
