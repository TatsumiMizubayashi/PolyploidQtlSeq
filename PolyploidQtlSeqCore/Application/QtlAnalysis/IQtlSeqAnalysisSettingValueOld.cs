﻿using PolyploidQtlSeqCore.QtlAnalysis;

namespace PolyploidQtlSeqCore.Application.QtlAnalysis
{
    /// <summary>
    /// QTL-Seq解析オプション値 インターフェイス
    /// </summary>
    public interface IQtlSeqAnalysisSettingValueOld : IQtlAnalysisScenarioSettingValue
    {
        /// <summary>
        /// InputVCFファイルPathを取得する。
        /// </summary>
        string InputVcf { get; }

        /// <summary>
        /// パラメータファイルPathを取得する。
        /// </summary>
        string ParameterFile { get; }
    }
}