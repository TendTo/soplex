# $Id: check.awk,v 1.18 2003/01/15 17:26:03 bzfkocht Exp $
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*   File....: check.awk                                                     *
#*   Name....: SoPlex Check Report Generator                                 *
#*   Author..: Thorsten Koch                                                 *
#*   Copyright by Author, All rights reserved                                *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
function abs(x)
{
   if (x < 0)
      return -x;
   if (x > 0)
      return x;
   return 0.0;  # get rid of -0.0
}
function printviol(x)
{
   if (x < 1e-9)
      printf("         ");
   else
      printf(" %.2e", abs(x));
}
BEGIN {
    print "$Id: check.awk,v 1.18 2003/01/15 17:26:03 bzfkocht Exp $";
    print "";
    line = "--------------------------------------------------------------------------------------------------------------------------------\n";
    printf(line);
    printf("Name         Rows   Cols Type   Iter     Time Objective                    maxVCons sumVCons maxVBoun sumVBoun maxVRedc sumVRedc\n");
}
/=opt=/          { sol[$2] = $3; }
/=type=/         { type = $2; }
/IEXAMP22/       { file = $5; }
/IEXAMP24/       { rows = $4; cols = $6; }
/IEXAMP27/       { time = $5; } 
/IEXAMP28/       { iter = $4; }
/IEXAMP29/       { obj  = $5; }
/IEXAMP31/       { infeas = 1; }
/IEXAMP32/       { infeas = 1; } 
/IEXAMP33/       { timeout = 1; }
/IEXAMP07/       { cvm = $4; cvs = $5; if (cvm > cvmax[type]) cvmax[type] = cvm; cvsum[type] += cvs; }
/IEXAMP09/       { bvm = $4; bvs = $5; if (bvm > bvmax[type]) bvmax[type] = bvm; bvsum[type] += bvs; }
/IEXAMP11/       { rcm = $5; rcs = $6; if (rcm > rcmax[type]) rcmax[type] = rcm; rcsum[type] += rcs; }
/=start=/        {
   type = "";
   for(i = 2; i <= NF; i++)
      type = type substr($i, 2);
}
/ready/       {
    n = split(file, a, "/");
    split(a[n], b, ".");
    name = b[1];

    if (sol[name] == "")
       print name, "nicht gefunden";
    else
    {
        if (name == prevname)
            printf("%25s", "");
        else
        {
            printf(line);
	    printf("%-10s %6d %6d ", name, rows, cols);
        }
	printf("%-3s %7d %8.2f ", type, iter, time);

        if (infeas)
	    printf("%-14s", "infeasible");
	else if (timeout)
	    printf("%-14s", "timeout");
	else
	    printf("%+e ", obj);

	if (timeout)
	   printf("\n");
	else
	{
            if (!infeas && sol[name] != "infeasible")
            {
                abserr = abs(sol[name] - obj);

                if (abs(sol[name]) >= 1e-5)
                    relerr = abserr / abs(sol[name]) * 1000.0
                else
                    relerr = 0.0;

    	        if ((abserr < 1e-4) || (relerr < 0.01))
		{
		   printf("ok            ");
		   pass[type]++;
		   passes++;
		}
		else
		{
		   printf("error %.2e", abserr);
		   fail[type]++;
		   fails++;
		}
		printviol(cvm);
		printviol(cvs);
		printviol(bvm);
		printviol(bvs);
		printviol(rcm);
		printviol(rcs); 
		print "";
	    }
	    else
	    {
	       if (infeas == 1 && sol[name] == "infeasible")
	       {
		  printf("ok\n");
		  pass[type]++;
		  passes++;
	       }
	       else
	       {
		  if (infeas && sol[name] != "infeasible")
		     printf("error %e\n", abs(sol[name]));
		  else
		     printf("error infeasible\n");
		  
		  fail[type]++;
		  fails++;
	       }
	    }
	}
        sum[type] += time;
        cnt[type]++;
        counts++;
        times += time;
    }
    prevname = name;
    timeout  = 0;
    infeas   = 0;
    obj      = 0;
    iter     = 0;
    time     = 0;
    rows     = 0;
    cols     = 0;
}
END {
    print "";
    printf(line);
    printf("Alg            Cnt  Pass  Fail       Time                                  maxVCons sumVCons maxVBoun sumVBoun maxVRedc sumVRedc\n");
    printf(line);
    for(i in sum)
    {
        printf("%-12s %5d %5d %5d %10.2f                                 ",
	   i, cnt[i], pass[i], fail[i], sum[i]);
	printviol(cvmax[i]);
	printviol(cvsum[i]);
	printviol(bvmax[i]);
	printviol(bvsum[i]);
	printviol(rcmax[i]);
	printviol(rcsum[i]);
	print "";

	cvmaxsum += cvmax[i];
	cvsumsum += cvsum[i];
	bvmaxsum += bvmax[i];
	bvsumsum += bvsum[i];
	rcmaxsum += rcmax[i];
	rcsumsum += rcsum[i];
    }
    printf(line);
    printf("%-12s %5d %5d %5d %10.2f                                 ",
       "Sum", counts, passes, fails, times);

    printviol(cvmaxsum);
    printviol(cvsumsum);
    printviol(bvmaxsum);
    printviol(bvsumsum);
    printviol(rcmaxsum);
    printviol(rcsumsum);
    print "";
    printf(line);
}






