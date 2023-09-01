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
            QtlAnalysisScenarioSettings.AddLongNameKeyValuePair(toLongNameDictionary);

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
            AnalysisChrSettings = new AnalysisChrSettings(optionValues, longNameParameterDictionary, userOptionDictionary);
            BcfToolsVariantCallSettings = new BcfToolsVariantCallSettings(optionValues, longNameParameterDictionary, userOptionDictionary);
            SnpEffSettings = new SnpEffSettings(optionValues, longNameParameterDictionary, userOptionDictionary);
            QtlAnalysisScenarioSettings = new QtlAnalysisScenarioSettings(optionValues, longNameParameterDictionary, userOptionDictionary);
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
        /// 解析染色体設定を取得する。
        /// </summary>
        public AnalysisChrSettings AnalysisChrSettings { get; }

        /// <summary>
        /// bcftools変異検出設定を取得する。
        /// </summary>
        public BcfToolsVariantCallSettings BcfToolsVariantCallSettings { get; }

        /// <summary>
        /// SnpEff設定を取得する。
        /// </summary>
        public SnpEffSettings SnpEffSettings { get; }

        /// <summary>
        /// QTL解析シナリオ設定を取得する。
        /// </summary>
        public QtlAnalysisScenarioSettings QtlAnalysisScenarioSettings { get; }

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
                .Concat(AnalysisChrSettings.ToParameterFileLines())
                .Concat(BcfToolsVariantCallSettings.ToParameterFileLines())
                .Concat(SnpEffSettings.ToParameterFileLines())
                .Concat(QtlAnalysisScenarioSettings.ToParameterFileLines());
            foreach (var line in parameterLineQuery)
            {
                writer.WriteLine(line);
            }
        }
    }
}
