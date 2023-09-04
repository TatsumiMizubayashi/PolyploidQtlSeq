﻿using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeqCore.QtlAnalysis;

namespace PolyploidQtlSeqCore.Application.QtlAnalysis
{
    /// <summary>
    /// QTL-Seq解析
    /// </summary>
    public class QtlSeqAnalysis
    {
        private readonly QtlSeqAnalysisSettings _option;

        /// <summary>
        /// QTL-Seq解析インスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">オプションの値</param>
        /// <param name="options">CommandOptions</param>
        public QtlSeqAnalysis(IQtlSeqAnalysisSettingValueOld optionValue, IReadOnlyCollection<CommandOption> options)
        {
            _option = new QtlSeqAnalysisSettings(optionValue, options);
        }

        /// <summary>
        /// QTL-seq解析を実行する。
        /// </summary>
        /// <returns>終了コード</returns>
        [Obsolete("削除予定")]
        public int Run()
        {
            var scenario = new QtlAnalysisScenario(_option.QtlAnalysisScenarioSettings);
            var code = scenario.Run(_option.InputVcf);

            var outputDir = _option.QtlAnalysisScenarioSettings.OutputDir;
            var paramsFilePath = outputDir.CreateFilePath("qtl.parameter.txt");
            _option.SaveParameterFile(paramsFilePath);

            return code;
        }

        private readonly QtlAnalysisScenarioSettings _settings;

        /// <summary>
        /// QTL-Seq解析インスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL解析シナリオ設定値</param>
        public QtlSeqAnalysis(IQtlAnalysisScenarioSettingValue settingValue)
        {
            _settings = new QtlAnalysisScenarioSettings(settingValue);
        }
    }
}
