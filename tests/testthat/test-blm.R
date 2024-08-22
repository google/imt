# Copyright 2024 Google LLC

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     https://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

test_that("The Bayesian Linear Model Works", {
  set.seed(8297)
  N <- 500
  fake_data <- tibble::tibble(
    y_pre = rnorm(n = N, mean = 11739, sd = 100),
    x1 = rnorm(n = N, mean = -200, sd = 10),
    x2 = rnorm(n = N, mean = 300, sd = 100),
    x3 = sample(c("a", "b", "c"), size = N, replace = TRUE),
    t = sample(c(FALSE, TRUE), size = N, replace = TRUE),
    y = rnorm(n = N, mean = 0, sd = 1) # This is just a placeholder for now
  )
  blm_fake <- imt::blm$new(
    y = "y", x = c("y_pre", "x1", "x2", "x3"),
    treatment = "t", data = fake_data, eta_mean = 0,
    eta_sd = 1, generate_fake_data = 1
  )
  expect_s3_class(blm_fake, "blm")
})
