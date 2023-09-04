using PolyploidQtlSeqCore.QtlAnalysis;

namespace PolyploidQtlSeqCore.Application.QtlAnalysis
{
    /// <summary>
    /// QTL-Seq解析
    /// </summary>
    public class QtlSeqAnalysis
    {
        private readonly QtlAnalysisScenarioSettings _settings;

        /// <summary>
        /// QTL-Seq解析インスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL解析シナリオ設定値</param>
        public QtlSeqAnalysis(IQtlAnalysisScenarioSettingValue settingValue)
        {
            _settings = new QtlAnalysisScenarioSettings(settingValue);
        }

        /// <summary>
        /// QTL-Seq解析を実行する。
        /// </summary>
        /// <param name="inputVcfPath">入力VCFファイルPath</param>
        /// <returns></returns>
        public int Run(string inputVcfPath)
        {
            var inputVcf = new InputVcf(inputVcfPath);
            var scenario = new QtlAnalysisScenario(_settings);
            var code = scenario.Run(inputVcf);

            return code;
        }
    }
}
