# Row-wise set operations on Ranges objects

Row-wise set operations on Ranges objects

## Usage

``` r
x %union% y

x %intersect% y

x %setdiff% y

between(x, y)

span(x, y)
```

## Arguments

- x, y:

  Ranges objects

## Value

A Ranges object

## Details

Each of these functions acts on the rows between pairs of Ranges object.
The function `%union%()`. will return the entire range between two
ranges objects assuming there are no gaps, if you would like to force
gaps use `span()` instead. The function `%intersect%()` will create a
new ranges object with a hit column indicating whether or not the two
ranges intersect. The function `%setdiff%()`will return the ranges for
each row in x that are not in the corresponding row of y. The function
`between()` will return the gaps between two ranges.

## See also

\[IRanges::punion()\]\[IRanges::pintersect()\]\[IRanges::pgap()\]\[IRanges::psetdiff()\]

## Examples

``` r
x <- as_iranges(data.frame(start = 1:10, width = 5))
# stretch x by 3 on the right
y <- stretch(anchor_start(x), 3)
# take the rowwise union
x %union% y
#> IRanges object with 10 ranges and 0 metadata columns:
#>            start       end     width
#>        <integer> <integer> <integer>
#>    [1]         1         8         8
#>    [2]         2         9         8
#>    [3]         3        10         8
#>    [4]         4        11         8
#>    [5]         5        12         8
#>    [6]         6        13         8
#>    [7]         7        14         8
#>    [8]         8        15         8
#>    [9]         9        16         8
#>   [10]        10        17         8
# take the rowwise intersection
x %intersect% y
#> IRanges object with 10 ranges and 0 metadata columns:
#>            start       end     width
#>        <integer> <integer> <integer>
#>    [1]         1         5         5
#>    [2]         2         6         5
#>    [3]         3         7         5
#>    [4]         4         8         5
#>    [5]         5         9         5
#>    [6]         6        10         5
#>    [7]         7        11         5
#>    [8]         8        12         5
#>    [9]         9        13         5
#>   [10]        10        14         5
# asymetric difference
y %setdiff% x
#> IRanges object with 10 ranges and 0 metadata columns:
#>            start       end     width
#>        <integer> <integer> <integer>
#>    [1]         6         8         3
#>    [2]         7         9         3
#>    [3]         8        10         3
#>    [4]         9        11         3
#>    [5]        10        12         3
#>    [6]        11        13         3
#>    [7]        12        14         3
#>    [8]        13        15         3
#>    [9]        14        16         3
#>   [10]        15        17         3
x %setdiff% y
#> IRanges object with 10 ranges and 0 metadata columns:
#>            start       end     width
#>        <integer> <integer> <integer>
#>    [1]         1         0         0
#>    [2]         2         1         0
#>    [3]         3         2         0
#>    [4]         4         3         0
#>    [5]         5         4         0
#>    [6]         6         5         0
#>    [7]         7         6         0
#>    [8]         8         7         0
#>    [9]         9         8         0
#>   [10]        10         9         0
# if there are gaps between the rows of each range use span
y <- as_iranges(data.frame(start = c(20:15, 2:5),
width = c(10:15,1:4)))
# fill in the gaps and take the rowwise union
span(x,y)
#> IRanges object with 10 ranges and 0 metadata columns:
#>            start       end     width
#>        <integer> <integer> <integer>
#>    [1]         1        29        29
#>    [2]         2        29        28
#>    [3]         3        29        27
#>    [4]         4        29        26
#>    [5]         5        29        25
#>    [6]         6        29        24
#>    [7]         2        11        10
#>    [8]         3        12        10
#>    [9]         4        13        10
#>   [10]         5        14        10
# find the gaps
between(x,y)
#> IRanges object with 10 ranges and 0 metadata columns:
#>            start       end     width
#>        <integer> <integer> <integer>
#>    [1]         6        19        14
#>    [2]         7        18        12
#>    [3]         8        17        10
#>    [4]         9        16         8
#>    [5]        10        15         6
#>    [6]        11        14         4
#>    [7]         3         6         4
#>    [8]         5         7         3
#>    [9]         7         8         2
#>   [10]         9         9         1
```
