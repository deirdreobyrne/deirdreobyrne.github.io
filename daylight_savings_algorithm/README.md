# Calculating Daylight Savings Time changeovers

## DST rules

There are [a limited number of rules](https://en.wikipedia.org/wiki/Daylight_saving_time_by_country) (16 at the time of writing) for when daylight savings
time comes into and out of effect. With the exception of the rules for Iran, the rules can be summarised as `[on|Fri before] [1st|2nd|last] [Thu|Fri|Sat|Sun]
in [Mar|Apr|Oct|Nov] at [time]`. Note, however, that DST rules can change, so it is prudent for a programmer to extend the possibilities, within reason. It
is also prudent to account for the fact that not all DST changes are plus one hour - at present, Lord Howe Island has a 30 minute DST change. DST has also
historically been 2 hours, especially during wartime.

## Step 1 - an algorithm to convert a date to a time interval since 1970

This algorithm is based on chapter 7 of Jean Meeus’ excellent 1991 book “Astronomical Algorithms” published by Willmann-Bell of Richmond, Virginia, USA.

```
// y - the year (>1970)
// m - the month (0 .. 11)
// d - the day of the month (1 .. 31)
// returns the number of days elapsed since 1970 Jan 1
int getDayNumber(int y, int m, int d) {
  int days;
  if (m < 2) {
    y--;
    m+=12;
  }
  days = (int)(y/100);
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

## Step 2 - Algorithm to implement the rules determining when the change occurs

Most computer languages have an implementation of something like `Date()`, from which the current year can be determined (e.g. by calling `Date().getYear()` or
equivalent), and the current number of milliseconds since 1970 Jan 1 (e.g. by calling `Date().getMillis()` or equivalent).

First we define a data structure
```
typedef struct {
  int daysPrior; // if the DST rule is "Fri before last Sun", then this is 2, otherwise it is zero
  int dayNumber; // 0 for 1st, 1 for 2nd, and 4 for last
  int dayOfWeek; // 0 for Sunday
  int month; // 0 for January
  float time; // the time of day (in hours) that the changeover happens
  bool isStart; // true if this represents the start (as opposed to end) of DST
} DSTRule;
```


```
// year - the year for which we are determining the start/end of DST
// rule - the DSTRule for the change we are considering
// dstTimezone - the "normal" timezone of the user in hours (positive east)
// dstOffset - the number of hours the DST changeover involves
// Returns the number of milliseconds since 1970 that the change happens
float getDstChangeTime(int year, DSTRule rule, float dstTimezone, float dstOffset) {
  int result;
  if (rule.dayNumber == 4) { // last X of this month? Work backwards from 1st of next month.
    if (++rule.month > 11) {
      year++;
      rule.month-=12;
    }
  }
  result = getDayNumber(year, rule.month, 1);
  // 1970 Jan 1 was Thursday, so (result % 7) is 0 for Thursday, hence ((result + 4) % 7) is 0 for Sunday
  if (rule.dayNumber == 4) {
    result -= 7 - (7 - ((result + 4) % 7) + rule.dayOfWeek) % 7;
  } else {
    result += 7 * rule.dayNumber + (14 + rule.dayOfWeek - ((result + 4) % 7)) % 7;
  }
  result -= rule.daysPrior;
  return 1000 * ((result * 86400) + (rule.time - dstTimezone - (rule.isStart ? 0 : dstOffset)) * 3600);
}
```

## Step 3 - algorithm to actually determine which changeover is next

```
typedef struct {
  float changeoverTime; // GMT millis since 1970 that the change happens
  DSTRule change; // The change that is due to happen
} DSTChange;

DSTChange getNextDSTChange(Date date, DSTRule dstStart, DSTRule dstEnd, float dstTimezone, float dstOffset) {
  float start = getDstChangeTime(date.getYear(), dstStart, dstTimezone, dstOffset);
  float end = getDstChangeTime(date.getYear(), dstEnd, dstTimezone, dstOffset);
  DSTChange result;
  if (start <= date.getTime()) {
    if (end <= date.getTime()) {
      // Both changes have happened for this year
      if (start < end) {
        // The start of DST is earlier than the end, so next change is a start of DST
        result.changeoverTime = getDstChangeTime(date.getYear()+1, dstStart, dstTimezone, dstOffset);
        result.change = dstStart;
      } else {
        // The end of DST is earlier than the start, so the next change is an end of DST
        result.changeoverTime = getDstChangeTime(date.getYear()+1, dstEnd, dstTimezone, dstOffset);
        result.change = dstEnd;
      }
    } else {
      result.changeoverTime = end;
      result.change = dstEnd;
    }
  } else {
    result.changeoverTime = start;
    result.change = dstStart;
  }
  return result;
}
```

