# Calculating Daylight Savings Time changeovers

## DST rules

There are [a limited number of rules](https://en.wikipedia.org/wiki/Daylight_saving_time_by_country) (16 at the time of writing) for when daylight savings
time comes into and out of effect. With the exception of the rules for Iran, the rules can be summarised as `[on|Fri before] [1st|2nd|last] [Thu|Fri|Sat|Sun]
in [Mar|Apr|Oct|Nov]`. Note, however, that DST rules can change, so it is prudent for a programmer to extend the possibilities, within reason.

## Step 1 - convert a date to a time interval since 1970

This algorithm is based on chapter 7 of Jean Meeus’ excellent 1991 book “Astronomical Algorithms” published by Willmann-Bell of Richmond, Virginia, USA.

```
// y - the year (>1970)
// m - the month (0 .. 11)
// d - the day of the month (1 .. 31)
// returns the number of days elapsed since 1970 Jan 1
int dayNumber(int y, int m, int d) {
  int days;
  if (m < 2) {
    y--;
    m+=12;
  }
  days = floor(y/100);
  days = 365*y + (y>>2) - days + (days>>2) + 30*m + (int)((3*m+6)/5) + d - 719531;
  return days;
}
```
The algorithm works by first re-framing the start of the year as 1st March. This makes it easy to deal with leap years, as a leap year just "starts" one day
later than normal. Then we calculate the number of the century, and store the result temporarily in the `days` variable. Finally we calculate the actual number of
days elapsed.

The elements of that final formula are -
- `365*y + (y>>2)` - first step in calculating the number of days - assuming every 4th year is a leap year
- `- days + (days>>2)` - remembering that `days` here actually refers to the century number, this calculates the Gregorian correction to the calendar
- `30*m + (int)((3*m+6)/5)` - calculate the number of days elapsed to the start of the month. This formula neatly accounts for the varying lengths of
the months when `m` is 2 for March and 13 for February. Note that the `((3*m+6)/5)` needs to be calculated using integer arithmetic. Meeus actually gives
this part of the formula as `floor(30.6001*(m+2))`, which is safer if your computer insists on floating-point.
- `d` - add the day of the month
- ` - 719531` - a correction so that 1970 Jan 1 is zero

