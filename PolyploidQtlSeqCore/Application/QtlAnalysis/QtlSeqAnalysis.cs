﻿using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeqCore.QtlAnalysis;

namespace PolyploidQtlSeqCore.Application.QtlAnalysis
{
    /// <summary>
    /// QTL-Seq解析
    /// </summary>
    public class QtlSeqAnalysis
    {
        private readonly QtlAnalysisCommandOption _option;

        /// <summary>
        /// QTL-Seq解析インスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">オプションの値</param>
        /// <param name="options">CommandOptions</param>
        public QtlSeqAnalysis(IQtlAnalysisCommandOptions optionValue, IReadOnlyCollection<CommandOption> options)
        {
            _option = new QtlAnalysisCommandOption(optionValue, options);
        }

        /// <summary>
        /// QTL-seq解析を実行する。
        /// </summary>
        /// <returns>終了コード</returns>
        public int Run()
        {
            var scenario = new QtlAnalysisScenario(_option.QtlAnalysisScenarioOptions);
            var code = scenario.Run(_option.InputVcf);

            var outputDir = _option.QtlAnalysisScenarioOptions.OutputDir;
            var paramsFilePath = outputDir.CreateFilePath("qtl.parameter.txt");
            _option.SaveParameterFile(paramsFilePath);

            return code;
        }
    }
}
