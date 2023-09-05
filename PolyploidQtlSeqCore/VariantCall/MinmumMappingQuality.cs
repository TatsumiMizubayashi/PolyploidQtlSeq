using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// Mapping Quality最小値
    /// </summary>
    public class MinmumMappingQuality
    {
        /// <summary>
        /// min-MQの最小値
        /// </summary>
        private const int MINIMUM = 0;

        /// <summary>
        /// min-MQの最大値
        /// </summary>
        private const int MAXIMUM = 60;

        /// <summary>
        /// min-MQの規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const int DEFAULT = 40;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "q";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "minMQ";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Minimum mapping quality at variant detection in bcftools mpileup.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        [Obsolete("削除予定")]
        public const string VALIDATION_ERROR_MESSAGE = "The -q option must be an integer greater than or equal to 0 and less than or equal to 60.";

        /// <summary>
        /// Mapping Qualityの最小値を作成する。
        /// </summary>
        /// <param name="minMq">MQ最小値</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public MinmumMappingQuality(int minMq, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, minMq, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// Mapping Qualityの最小値を作成する。
        /// </summary>
        /// <param name="minMq">MQ最小値</param>
        public MinmumMappingQuality(int minMq)
        {
            if (minMq < MINIMUM || minMq > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(minMq));

            Value = minMq;
        }

        /// <summary>
        /// MQ最小値を取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        internal string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Value}";
        }

    }
}
