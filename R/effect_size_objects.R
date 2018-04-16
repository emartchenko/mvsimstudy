#' An S4 class to represent the effect size the user specifies is present
#' between the treatment group and the control group in the scenario.
#'
#' @slot effects The vector of differences in mean between the two groups
#' being considered
#' @slot id Numeric assigned by the check_effects method when the object is
#' created. Unique for each effects vector created.
#' @export
setClass(Class = "effect_size", representation(effects = "vector",
    id = "numeric"))

# A separate environment for storing list of effect size objects and returning
# the index for a previously declared effect_size if the user attempts to
# created an identical object
effect_size.env <- new.env()

#' @title Check effect size vector
#' @description Checks the vector specifying the effects used for the scenario
#' to ensure no duplicate effect sizes are being stored.
#' @param es_vector A vector of user-specified effect sizes.
#' @return If vector has already been added previously the index for
#' the effect_size object location in the list is return. If it has not, then
#' the vector is added to the list and 0 returned to indicate that a new
#' effect_size object should be created.
#' @details An internal function that carries out it's operation in the
#' environment created for effect_size checking. Keeps a list of vectors that
#' were specified by the user, as well as a list of effect_size objects in this
#' environment.
#' @rdname check_effects
#'
check_effects <- local(function(es_vector) {

    if (!exists("list_sizes")) {
        list_sizes <<- list()
    }

    where <- Position(function(x) identical(x, es_vector), list_sizes,
        nomatch = 0)

    if (where == 0) {
        list_sizes <<- c(list_sizes, list(es_vector))
    }

    return(where)
}, env = effect_size.env)

#' @title effect_size object constructor
#' @description Creates (or if it already exists, simply returns) a unique effect size object
#' @param effects A vector of the effect size the user would like to apply to
#' the simulation they are running
#' @return An effect size object that can be passed to the simulation functions
#' @details Each effect size object is unique. Duplicates will not be created.
#' @examples
#' effect_size(c(rep(0,5), 0.4, rep(0.3, 6)))
#' effect_size(c(0.3, 0.4, 0, 0, 0, 0.2))
#' @rdname effect_size
#' @export
effect_size <- local(function(effects) {
    # check if this effect size has already been created
    result <- check_effects(effects)
    # create a list of effect_sizes if it doesn't already exist
    if (!exists("list_effect_objects")) {
        list_effect_objects <<- list()
    }
    # if not duplicate, create a new effect_size object with this
    # vector
    if (result == 0) {
        length_list <- length(list_effect_objects)
        obj <- new("effect_size", effects = effects, id = length_list +
            1)
        list_effect_objects <<- c(list_effect_objects, obj)
    } else {
        # if duplicate return the object with this effect size
        obj <- list_effect_objects[[result]]
    }

    return(obj)
}, envir = effect_size.env)

