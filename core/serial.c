static int treebase( int *nlen, char **aseq, int nadd, char *mergeoralign, char **mseq1, char **mseq2, int ***topol, Treedep *dep, int **memhist, double ***cpmxhist, double *effarr, double **newdistmtx, int *selfscore, int *alloclen, int (*callback)(int, int, char*) )
{
	int l, len1, len2, i, m, immin, immax;
	int len1nocommongap, len2nocommongap;
	int clus1, clus2;
	double pscore, tscore;
	char *indication1 = NULL, *indication2 = NULL;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	int *fftlog = NULL; // fixed at 2006/07/26
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int ffttry;
	int m1, m2;
	int *gaplen = NULL;
	int *gapmap = NULL;
	int *alreadyaligned = NULL;
	double **dynamicmtx = NULL;
	double ssi, ssm, bunbo;
	int tm, ti;
	int gapmaplen;
	int **localmem = NULL;
	double **cpmxchild0, **cpmxchild1;
	double orieff1, orieff2;
#if SKIP
	int **skiptable1 = NULL, **skiptable2 = NULL;
#endif
#if 0
	int i, j;
#endif


	if( effarr1 == NULL ) 
	{
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		fftlog = AllocateIntVec( njob );
		gaplen = AllocateIntVec( *alloclen+10 );
		gapmap = AllocateIntVec( *alloclen+10 );
		alreadyaligned = AllocateIntVec( njob );
		if( specificityconsideration )
			dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );
		localmem = calloc( sizeof( int * ), 2 );
	}
	for( i=0; i<njob-nadd; i++ ) alreadyaligned[i] = 1;
	for( i=njob-nadd; i<njob; i++ ) alreadyaligned[i] = 0;

	if( callback && callback( 0, 50, "Progressive alignment" ) ) goto chudan_tbfast;

	for( l=0; l<njob; l++ ) fftlog[l] = 1;

#if 0 // chain you
	localmem[0][0] = -1;
	localmem[1][0] = -1;
	clus1 = 1;// chain ni hitsuyou
#endif

#if 0
	reporterr(       "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		reporterr(       "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

//	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );


//	writePre( njob, name, nlen, aseq, 0 );

	tscore = 0.0;
	for( l=0; l<njob-1; l++ )
	{
		m1 = topol[l][0][0];
		m2 = topol[l][1][0];
//		reporterr( " at the beginning of the loop, clus1,clus2=%d,%d\n", clus1, clus2 );

//		reporterr( "l=%d, dep[l].child0=%d, dep[l].child1=%d\n", l, dep[l].child0, dep[l].child1 );
		if( dep[l].child0 == -1 ) cpmxchild0 = NULL; else cpmxchild0 = cpmxhist[dep[l].child0];
		if( dep[l].child1 == -1 ) cpmxchild1 = NULL; else cpmxchild1 = cpmxhist[dep[l].child1];
//		reporterr( "cpmxchild0=%p, cpmxchild1=%p\n", cpmxchild0, cpmxchild1 );

#if 0
		if(  l > 0 && dep[l].child0 == l-1 && dep[l].child1 == -1 && dep[dep[l].child0].child1 == -1 )
		{
			localmem[0][clus1] = topol[l-1][1][0];
			localmem[0][clus1+1] = -1;

			localmem[1][0] = topol[l][1][0];
			localmem[1][1] = -1;
		}
		else
		{
			localmem[0][0] = -1;
			posinmem = topolorderz( localmem[0], topol, dep, l, 0 ) - localmem[0];
			localmem[1][0] = -1;
			posinmem = topolorderz( localmem[1], topol, dep, l, 1 ) - localmem[1];
		}
#else
		if( dep[l].child0 == -1 ) 
		{
			localmem[0] = calloc( sizeof( int ), 2 );
			localmem[0][0] = m1;
			localmem[0][1] = -1;
			clus1 = 1;
		}
		else
		{
			localmem[0] = memhist[dep[l].child0];
			clus1 = intlen( localmem[0] );
		}
		if( dep[l].child1 == -1 ) 
		{
			localmem[1] = calloc( sizeof( int ), 2 );
			localmem[1][0] = m2;
			localmem[1][1] = -1;
			clus2 = 1;
		}
		else
		{
			localmem[1] = memhist[dep[l].child1];
			clus2 = intlen( localmem[1] );
		}

		if( l != njob-2 )
		{
			memhist[l] = calloc( sizeof( int ), clus1+clus2+1 );
			intcpy( memhist[l], localmem[0] );
			intcpy( memhist[l]+clus1, localmem[1] );
			memhist[l][clus1+clus2] = -1;
		}
#endif

		if( mergeoralign[l] == 'n' )
		{
//			reporterr(       "SKIP!\n" );
//			free( topol[l][0] ); topol[l][0] = NULL;
//			free( topol[l][1] ); topol[l][1] = NULL;
//			free( topol[l] ); topol[l] = NULL;
			continue;
		}

//		reporterr(       "\ndistfromtip = %f\n", dep[l].distfromtip );
		if( specificityconsideration )
			makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );
		else
			dynamicmtx = n_dis_consweight_multi;
//		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, ( dep[l].distfromtip - 0.2 ) * 3 );


		len1 = strlen( aseq[m1] );
		len2 = strlen( aseq[m2] );
		if( *alloclen < len1 + len2 )
		{
			reporterr(       "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  ); 
			gaplen = realloc( gaplen, ( *alloclen + 10 ) * sizeof( int ) );
			if( gaplen == NULL )
			{
				reporterr(       "Cannot realloc gaplen\n" );
				exit( 1 );
			}
			gapmap = realloc( gapmap, ( *alloclen + 10 ) * sizeof( int ) );
			if( gapmap == NULL )
			{
				reporterr(       "Cannot realloc gapmap\n" );
				exit( 1 );
			}
			reporterr(       "done. *alloclen = %d\n", *alloclen );
		}

#if 1 // CHUUI@@@@
		clus1 = fastconjuction_noname( localmem[0], aseq, mseq1, effarr1, effarr, indication1, 0.0, &orieff1 );
		clus2 = fastconjuction_noname( localmem[1], aseq, mseq2, effarr2, effarr, indication2, 0.0, &orieff2 );
#else
		clus1 = fastconjuction_noname( topol[l][0], aseq, mseq1, effarr1, effarr, indication1, 0.0 );
		clus2 = fastconjuction_noname( topol[l][1], aseq, mseq2, effarr2, effarr, indication2, 0.0 );
//		clus1 = fastconjuction_noweight( topol[l][0], aseq, mseq1, effarr1,  indication1 );
//		clus2 = fastconjuction_noweight( topol[l][1], aseq, mseq2, effarr2,  indication2 );
#endif











		if( mergeoralign[l] == '1' || mergeoralign[l] == '2' )
		{
			newgapstr = "=";
		}
		else
			newgapstr = "-";

		len1nocommongap = len1;
		len2nocommongap = len2;
		if( mergeoralign[l] == '1' ) // nai
		{
			findcommongaps( clus2, mseq2, gapmap );
			commongappick( clus2, mseq2 );
			len2nocommongap = strlen( mseq2[0] );
		}
		else if( mergeoralign[l] == '2' )
		{
			findcommongaps( clus1, mseq1, gapmap );
			commongappick( clus1, mseq1 );
			len1nocommongap = strlen( mseq1[0] );
		}

#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
#endif


#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
#endif


//		fprintf( trap_g, "\nSTEP-%d\n", l );
//		fprintf( trap_g, "group1 = %s\n", indication1 );
//		fprintf( trap_g, "group2 = %s\n", indication2 );

//		reporterr(       "\rSTEP % 5d / %d %d-%d", l+1, njob-1, clus1, clus2 );
		if( l < 500 || l % 100 == 0 ) reporterr(       "\rSTEP % 5d / %d ", l+1, njob-1 );
		if( callback && callback( 0, 50+50*l/(njob-1), "Progressive alignment" ) ) goto chudan_tbfast;
#if 0
		reporterr( "\nclus1=%d, clus2=%d\n", clus1, clus2 );
#endif

#if 0
		reporterr(       "STEP %d /%d\n", l+1, njob-1 );
		reporterr(       "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
		reporterr(       "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
#endif

/*
		reporterr(       "before align all\n" );
		display( aseq, njob );
		reporterr(       "\n" );
		reporterr(       "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		reporterr(       "\n" );
		reporterr(       "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		reporterr(       "\n" );
*/


		if( !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
		{
			reporterr(       "\nlen1=%d, len2=%d, Switching to the memsave mode\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}

//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
		else						   ffttry = 0;
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000); // v6.708
//		reporterr(       "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (double)len1/fftlog[m1], clus1, (double)len2/fftlog[m2], clus2 );

		if( force_fft || ( use_fft && ffttry ) )
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " f\b\b" );
			if( alg == 'M' )
			{
				if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
			{
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
//				reporterr(       "######### mseq1[0] = %s\n", mseq1[0] );
			}
		}
		else
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " d\b\b" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
					if( l < 500 || l % 100 == 0 ) if( cpmxchild1 || cpmxchild0 ) reporterr(       " h" );
//					reporterr(       "%d-%d", clus1, clus2 );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					break;
				case( 'd' ):
					if( 1 && clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = D__align_ls( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					}
					break;
				case( 'A' ):
					if( clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
						if( l < 500 || l % 100 == 0 ) if( cpmxchild1 || cpmxchild0 ) reporterr(       " h" );
//						reporterr(       "\n\n %d - %d (%d x %d) : \n", topol[l][0][0], topol[l][1][0], clus1, clus2 );
						pscore = A__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, localmem[0][0], 1, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					}

					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}
#if SCOREOUT
		reporterr(       "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

//		writePre( njob, name, nlen, aseq, 0 );

		if( disp ) display( aseq, njob );
//		reporterr(       "\n" );

		if( mergeoralign[l] == '1' ) // jissainiha nai. atarashii hairetsu ha saigo dakara.
		{
			reporterr( "Check source!!!\n" );
			exit( 1 );
		}
		if( mergeoralign[l] == '2' )
		{
//			if( localkeeplength ) ndeleted += deletenewinsertions( clus1, clus2, mseq1, mseq2, NULL );
//			for( i=0; i<clus1; i++ ) reporterr(       ">STEP0 mseq1[%d] = \n%s\n", i, mseq1[i] );
//			for( i=0; i<clus2; i++ ) reporterr(       ">STEP0 mseq2[%d] = \n%s\n", i, mseq2[i] );
			gapmaplen = strlen( mseq1[0] )-len1nocommongap+len1;
			adjustgapmap( gapmaplen, gapmap, mseq1[0] );
#if 0
			reporterr( "\n" );
			for( i=0; i<clus1; i++ ) reporterr(       ">STEP1 mseq1[%d] = \n%s\n", i, mseq1[i] );
			for( i=0; i<clus2; i++ ) reporterr(       ">STEP1 mseq2[%d] = \n%s\n", i, mseq2[i] );
#endif
//			if( clus1 + clus2 < njob ) restorecommongaps( njob, aseq, topol[l][0], topol[l][1], gapmap, *alloclen, '-' );
			if( smoothing )
			{
				restorecommongapssmoothly( njob, njob-(clus1+clus2), aseq, localmem[0], localmem[1], gapmap, *alloclen, '-' );
				findnewgaps( clus1, 0, mseq1, gaplen );
				insertnewgaps_bothorders( njob, alreadyaligned, aseq, localmem[0], localmem[1], gaplen, gapmap, gapmaplen, *alloclen, alg, '-' );
			}
			else
			{
				restorecommongaps( njob, njob-(clus1+clus2), aseq, localmem[0], localmem[1], gapmap, *alloclen, '-' );
				findnewgaps( clus1, 0, mseq1, gaplen );
				insertnewgaps( njob, alreadyaligned, aseq, localmem[0], localmem[1], gaplen, gapmap, *alloclen, alg, '-' );
			}

#if 0
			reporterr( "\n" );
			for( i=0; i<clus1; i++ ) reporterr(       ">STEP3 mseq1[%d] = \n%s\n", i, mseq1[i] );
			for( i=0; i<clus2; i++ ) reporterr(       ">STEP3 mseq2[%d] = \n%s\n", i, mseq2[i] );
#endif

#if 0
			for( i=0; i<njob; i++ ) eq2dash( aseq[i] );
			for( i=0; i<clus1; i++ ) 
			{
				reporterr( "mseq1[%d] bef change = %s\n", i, mseq1[i] );
				eq2dash( mseq1[i] );
				reporterr( "mseq1[%d] aft change = %s\n", i, mseq1[i] );
			}
			for( i=0; i<clus2; i++ ) 
			{
				reporterr( "mseq2[%d] bef change = %s\n", i, mseq2[i] );
				eq2dash( mseq2[i] );
				reporterr( "mseq2[%d] aft change = %s\n", i, mseq2[i] );
			}
			for( i=0; i<clus1; i++ ) eq2dash( mseq1[i] );
			for( i=0; i<clus2; i++ ) eq2dash( mseq2[i] );
#endif


			eq2dashmatometehayaku( mseq1, clus1 );
			eq2dashmatometehayaku( mseq2, clus2 );

			for( i=0; (m=localmem[1][i])>-1; i++ ) alreadyaligned[m] = 1;
		}

		if( newdistmtx ) // tsukawanai
		{
#if 0
			reporterr( "group1 = " );
			for( i=0; i<clus1; i++ ) reporterr( "%d ", topol[l][0][i] );
			reporterr( "\n" );
			reporterr( "group2 = " );
			for( m=0; m<clus2; m++ ) reporterr( "%d ", topol[l][1][m] );
			reporterr( "\n" );
#endif
#if SKIP
			skiptable1 = AllocateIntMtx( clus1, 0 );
			skiptable2 = AllocateIntMtx( clus2, 0 );
			makeskiptable( clus1, skiptable1, mseq1 ); // allocate suru.
			makeskiptable( clus2, skiptable2, mseq2 ); // allocate suru.
#endif
			for( i=0; i<clus1; i++ ) 
			{
#if SKIP
//				makeskiptable( 1, skiptable1, mseq1+i ); // allocate suru.
#endif
				ti = localmem[0][i];
				ssi = selfscore[localmem[0][i]];
				for( m=0; m<clus2; m++ )
				{
					ssm = selfscore[localmem[1][m]];
					tm = localmem[1][m];
					if( ti<tm )
					{
						immin = ti;
						immax = tm;
					}
					else
					{
						immin = tm;
						immax = ti;
					}
					bunbo = MIN( ssi, ssm );
					if( bunbo == 0.0 )
						newdistmtx[immin][immax-immin] = 2.0; // 2013/Oct/17
					else
#if SKIP
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscorefast( mseq1[i], mseq2[m], skiptable1[i], skiptable2[m], penalty_dist ) / bunbo ) * 2.0;
#else
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscore11( mseq1[i], mseq2[m], penalty_dist ) / bunbo ) * 2.0;
#endif
				}
			}
#if SKIP
			FreeIntMtx( skiptable1 ); skiptable1 = NULL;
			FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif
		}

//		free( topol[l][0] ); topol[l][0] = NULL;
//		free( topol[l][1] ); topol[l][1] = NULL;
//		free( topol[l] ); topol[l] = NULL;


//		reporterr(       ">514\n%s\n", aseq[514] );
		free( localmem[0] );
		free( localmem[1] );
	}

#if SCOREOUT
	reporterr(       "totalscore = %10.2f\n\n", tscore );
#endif
	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
	D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	free( effarr1 );
	free( effarr2 );
	free( indication1 );
	free( indication2 );
	free( fftlog );
	free( gaplen );
	free( gapmap );
	if( specificityconsideration )
		FreeDoubleMtx( dynamicmtx );
	free( alreadyaligned );
	free( localmem );
	effarr1 = NULL;
	return( 0 );

	chudan_tbfast:

	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
	D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	if( effarr1 ) free( effarr1 ); effarr1 = NULL;
	if( effarr2 ) free( effarr2 ); effarr2 = NULL;
	if( indication1 ) free( indication1 ); indication1 = NULL;
	if( indication2 ) free( indication2 ); indication2 = NULL;
	if( fftlog ) free( fftlog ); fftlog = NULL;
	if( gaplen ) free( gaplen ); gaplen = NULL;
	if( gapmap ) free( gapmap ); gapmap = NULL;
	if( alreadyaligned ) free( alreadyaligned ); alreadyaligned = NULL;
	if( specificityconsideration )
	{
		if( dynamicmtx ) FreeDoubleMtx( dynamicmtx ); dynamicmtx = NULL;
	}
	if( localmem ) free( localmem ); localmem = NULL;
#if SKIP
	if( skiptable1 ) FreeIntMtx( skiptable1 ); skiptable1 = NULL;
	if( skiptable2 ) FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif

	return( 1 );
}
