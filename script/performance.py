import numpy as np

"""
Author: Fabiana Patalano
Mail: fabiana.patalano97@gmail.com
Last updated: 11/05/2021
"""


def relaxation(true_value, prediction, n=0):
    """

    :param true_value:
    :param prediction:
    :param n: relax parameter
    :return:
    """

    if np.isnan(true_value).all():
        raise Exception('ValueError')

    TP = 0
    TN = 0
    FP = 0
    FN = 0
    tot= []

    padded_true_interd = np.pad(true_value, ((n, n), (n, n)), 'constant', constant_values=np.nan)
    predicted_contact = np.where(~np.isnan(prediction) & (prediction != 20))
    predicted_non_contact = np.where(np.isnan(prediction) | (prediction == 20))
    for pred_cont_coord in zip(predicted_contact[0], predicted_contact[1]):
        if not np.isnan(padded_true_interd[pred_cont_coord[0]:pred_cont_coord[0] + (2 * n + 1),
                        pred_cont_coord[1]:pred_cont_coord[1] + (2 * n + 1)]).all():
            TP += 1
            tot.append('TP')
        elif np.isnan(padded_true_interd[pred_cont_coord]):
            FP += 1
            tot.append('FP')
    for pred_non_cont_coord in zip(predicted_non_contact[0], predicted_non_contact[1]):
        if np.isnan(padded_true_interd[pred_non_cont_coord]):
            TN += 1
            tot.append('TN')
        elif not np.isnan(padded_true_interd[pred_non_cont_coord]):
            FN += 1
            tot.append('FN')
    ppv = TP / (TP + FP)
    return TP, FP, TN, FN,ppv, tot
