using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeqCore.Options;
using PolyploidQtlSeqCore.QtlAnalysis;

namespace PolyploidQtlSeqCore.Application.QtlAnalysis
{
    /// <summary>
    /// QtlSeq解析設定
    /// </summary>
    [Obsolete("オプションスイッチ周りの機能を削除する予定")]
    internal class QtlSeqAnalysisSettings
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary;

        /// <summary>
        /// static コンストラクタ
        /// </summary>
        static QtlSeqAnalysisSettings()
        {
            var toLongNameDictionary = new Dictionary<string, string>()
            {
                [InputVcf.SHORT_NAME] = InputVcf.LONG_NAME,
                [InputVcf.LONG_NAME] = InputVcf.LONG_NAME
            };

            QtlAnalysisScenarioOptions.AddLongNameKeyValuePair(toLongNameDictionary);
            _toLongNameDictionary = toLongNameDictionary;
        }
                

        /// <summary>
        /// QTL解析コマンドオプションを作成する。
        /// </summary>
        /// <param name="optionValues">オプションの値</param>
        /// <param name="options">CommandOptions</param>
        public QtlSeqAnalysisSettings(IQtlSeqAnalysisSettingValue optionValues, IReadOnlyCollection<CommandOption> options)
        {
            ParameterFile = new ParameterFileParser(optionValues.ParameterFile);
            var longNameParameterDictionary = ParameterFile.ToParameterDictionary(_toLongNameDictionary);
            var userOptionDictionary = UserSpecifiedLongNameDictionaryCreator.Create(options);

            InputVcf = new InputVcf(optionValues.InputVcf, longNameParameterDictionary, userOptionDictionary);
            QtlAnalysisScenarioOptions = new QtlAnalysisScenarioOptions(optionValues, longNameParameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// 入力VCFファイル
        /// </summary>
        public InputVcf InputVcf { get; }

        /// <summary>
        /// QTL解析シナリオオプション
        /// </summary>
        public QtlAnalysisScenarioOptions QtlAnalysisScenarioOptions { get; }

        /// <summary>
        /// パラメータファイル
        /// </summary>
        public ParameterFileParser ParameterFile { get; }

        /// <summary>
        /// パラメータファイルを保存する。
        /// </summary>
        /// <param name="filePath">パラメータファイルPath</param>
        public void SaveParameterFile(string filePath)
        {
            using var writer = new StreamWriter(filePath);

            writer.WriteLine("#qtl Command");
            writer.WriteLine("#LongName\tValue");
            writer.WriteLine(InputVcf.ToParameterFileLine());
            
            foreach(var line in QtlAnalysisScenarioOptions.ToParameterFileLines())
            {
                writer.WriteLine(line);
            }
        }
    }
}
