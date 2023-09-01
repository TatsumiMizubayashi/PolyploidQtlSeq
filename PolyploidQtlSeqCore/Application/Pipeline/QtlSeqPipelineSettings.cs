using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeqCore.Mapping;
using PolyploidQtlSeqCore.Options;
using PolyploidQtlSeqCore.QtlAnalysis;
using PolyploidQtlSeqCore.QtlAnalysis.Chr;
using PolyploidQtlSeqCore.Share;
using PolyploidQtlSeqCore.VariantCall;

namespace PolyploidQtlSeqCore.Application.Pipeline
{
    /// <summary>
    /// QtlSeqパイプライン設定
    /// </summary>
    internal class QtlSeqPipelineSettings
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary;

        /// <summary>
        /// static コンストラクタ
        /// </summary>
        static QtlSeqPipelineSettings()
        {
            var toLongNameDictionary = new Dictionary<string, string>()
            {
                [ReferenceSequence.SHORT_NAME] = ReferenceSequence.LONG_NAME,
                [ReferenceSequence.LONG_NAME] = ReferenceSequence.LONG_NAME,
            };

            MappingSampleSettings.AddLongNameKeyValuePair(toLongNameDictionary);
            AnalysisChrSettings.AddLongNameKeyValuePair(toLongNameDictionary);
            BcfToolsVariantCallSettings.AddLongNameKeyValuePair(toLongNameDictionary);
            SnpEffSettings.AddLongNameKeyValuePair(toLongNameDictionary);
            QtlAnalysisScenarioOptions.AddLongNameKeyValuePair(toLongNameDictionary);

            _toLongNameDictionary = toLongNameDictionary;
        }

        /// <summary>
        /// QtlSeqパイプライン設定を作成する。
        /// </summary>
        /// <param name="optionValues">オプションの値</param>
        /// <param name="options">コマンドオプション</param>
        public QtlSeqPipelineSettings(IQtlSeqPipelineSettingValue optionValues, IReadOnlyCollection<CommandOption> options)
        {
            ParameterFile = new ParameterFileParser(optionValues.ParameterFile);
            var longNameParameterDictionary = ParameterFile.ToParameterDictionary(_toLongNameDictionary);
            var userOptionDictionary = UserSpecifiedLongNameDictionaryCreator.Create(options);

            //ReferenceSequence = new ReferenceSequence(optionValues.ReferenceSequence, longNameParameterDictionary, userOptionDictionary);

            MappingSettings = new MappingSettings(optionValues);
            MappingSampleSettings = new MappingSampleSettings(optionValues, longNameParameterDictionary, userOptionDictionary);
            AnalysisChrOption = new AnalysisChrSettings(optionValues, longNameParameterDictionary, userOptionDictionary);
            BcfToolsVariantCallOption = new BcfToolsVariantCallSettings(optionValues, longNameParameterDictionary, userOptionDictionary);
            SnpEffOption = new SnpEffSettings(optionValues, longNameParameterDictionary, userOptionDictionary);
            QtlAnalysisScenarioOptions = new QtlAnalysisScenarioOptions(optionValues, longNameParameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// リファレンスシークエンスを取得する。
        /// </summary>
        [Obsolete("削除予定")]
        public ReferenceSequence ReferenceSequence { get; }

        /// <summary>
        /// Mapping設定を取得する。
        /// </summary>
        public MappingSettings MappingSettings { get; }

        /// <summary>
        /// Mappingサンプル設定を取得する。
        /// </summary>
        public MappingSampleSettings MappingSampleSettings { get; }

        /// <summary>
        /// 解析染色体オプションを取得する。
        /// </summary>
        public AnalysisChrSettings AnalysisChrOption { get; }

        /// <summary>
        /// bcftools変異検出オプションを取得する。
        /// </summary>
        public BcfToolsVariantCallSettings BcfToolsVariantCallOption { get; }

        /// <summary>
        /// SnpEffオプションを取得する。
        /// </summary>
        public SnpEffSettings SnpEffOption { get; }

        /// <summary>
        /// QTL解析シナリオオプションを取得する。
        /// </summary>
        public QtlAnalysisScenarioOptions QtlAnalysisScenarioOptions { get; }

        /// <summary>
        /// パラメータファイルを取得する。
        /// </summary>
        public ParameterFileParser ParameterFile { get; }

        /// <summary>
        /// パラメータファイルを保存する。
        /// </summary>
        /// <param name="filePath">パラメータファイルPath</param>
        [Obsolete("削除予定")]
        public void SaveParameterFile(string filePath)
        {
            using var writer = new StreamWriter(filePath);

            writer.WriteLine("#QTL-Seq Command");
            writer.WriteLine("#LongName\tValue");
            //writer.WriteLine(ReferenceSequence.ToParameterFileLine());

            var parameterLineQuery = MappingSampleSettings.ToParameterFileLines()
                .Concat(AnalysisChrOption.ToParameterFileLines())
                .Concat(BcfToolsVariantCallOption.ToParameterFileLines())
                .Concat(SnpEffOption.ToParameterFileLines())
                .Concat(QtlAnalysisScenarioOptions.ToParameterFileLines());
            foreach (var line in parameterLineQuery)
            {
                writer.WriteLine(line);
            }
        }
    }
}
