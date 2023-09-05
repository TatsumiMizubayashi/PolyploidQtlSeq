using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// BaseQuality最低値
    /// </summary>
    public class MinmumBaseQuality
    {
        /// <summary>
        /// min-BQの最小値
        /// </summary>
        public const int MINIMUM = 0;

        /// <summary>
        /// min-BQの最大値
        /// </summary>
        public const int MAXIMUM = 60;

        /// <summary>
        /// min-BQの規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const int DEFAULT = 13;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "Q";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "minBQ";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Minimum base quality at variant detection in bcftools mpileup.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        [Obsolete("削除予定")]
        public const string VALIDATION_ERROR_MESSAGE = "The -Q option must be an integer greater than or equal to 0 and less than or equal to 60.";

        /// <summary>
        /// Base Qualityの最小値を作成する。
        /// </summary>
        /// <param name="minBq">Base Quality最小値</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public MinmumBaseQuality(int minBq, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, minBq, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// Base Qualityの最小値を作成する。
        /// </summary>
        /// <param name="minBq">Base Quality最小値</param>
        public MinmumBaseQuality(int minBq)
        {
            if (minBq < MINIMUM || minBq > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(minBq));

            Value = minBq;
        }

        /// <summary>
        /// Base Quality最小値を取得する。
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
