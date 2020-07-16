/* VHS Pawprints or Tiles  */
SELECT 
  *
FROM 
  vistaqc 
WHERE 
  is_tile = 'True' 
  /* is_stack = 'True' */
  /* VHS */
  AND prog = '179.A-2010'
  /* AND prog like '091.A-0426' */
  /* AND (nightobs between 20120401 and 20120930) */
  /* AND (nightobs between 20131001 and 20130331) */
  /* AND Semester = 'B' */
  /* AND (nightobs between 20100120 and 20100930) */
