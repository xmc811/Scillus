
context("Get colors")
library(Scillus)

test_that("Correct output of get_colors", {
        expect_equal(get_colors(v = c(1,3), pal = "Set1"), c("#E41A1C", "#4DAF4A"))
        expect_error(get_colors(v = -1))
})
